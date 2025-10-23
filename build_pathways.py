#!/usr/bin/env python3
"""
BioCyc to WikiPathways GPML Converter
======================================

Script for converting BioCyc pathway data to GPML format.

Usage:
    python build_pathways.py <data_dir> <output_dir> [options]

Arguments:
    data_dir            : Directory containing BioCyc .dat files (compounds.dat, genes.dat, etc.)
    output_dir          : Directory where GPML pathway files will be saved

Options:
    --include-reactions : Also build single reaction files for unused reactions
    --layout grid       : Use grid layout (default)
    --layout forceatlas2: Use ForceAtlas2 force-directed layout

Examples:
    python build_pathways.py ./data ./output
    python build_pathways.py ./data ./output --include-reactions
    python build_pathways.py ./data ./output --layout forceatlas2
    python build_pathways.py ./data ./output --include-reactions --layout forceatlas2
"""

import sys
import io
import os
import re

# Fix UTF-8 encoding issues on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from scripts.build_functions.general_pathwaybuilder import build_individual_pathways
from scripts.build_functions.pathway_builder_core import CompletePathwayBuilderWithGenes
from scripts.data_structure.wiki_data_structure import Pathway, Graphics, Xref, Author, Comment, Property as PathwayProperty
from scripts.parsing_functions import parsing_utils
from datetime import datetime


def find_unused_reactions(builder, reactions_file, pathways_file):
    """Find all reactions that are NOT used in any pathway."""
    # Get all reactions from reactions.dat
    reactions_processor = parsing_utils.read_and_parse(reactions_file)
    all_reactions = set()

    for record in reactions_processor.records:
        reaction_id = record.get('UNIQUE-ID')
        if reaction_id:
            all_reactions.add(reaction_id)

    # Get reactions used in pathways
    pathways_processor = parsing_utils.read_and_parse(pathways_file)
    reactions_in_pathways = set()

    for record in pathways_processor.records:
        reaction_list = record.get('REACTION-LIST', [])
        if not isinstance(reaction_list, list):
            reaction_list = [reaction_list] if reaction_list else []
        reactions_in_pathways.update(reaction_list)

    # Find unused reactions
    unused_reactions = all_reactions - reactions_in_pathways
    return sorted(list(unused_reactions))


def get_reaction_info(builder, reaction_id):
    """Get reaction information from the reaction record."""
    for record in builder.reaction_processor.records:
        if record.get('UNIQUE-ID') == reaction_id:
            return record
    return None


def build_single_reaction_pathway(builder, reaction_id):
    """Build a pathway object containing a single reaction."""
    # Get reaction information
    reaction_record = get_reaction_info(builder, reaction_id)

    # Get reaction name and EC number
    reaction_name = reaction_id
    ec_number = None
    organism = "TAX-3702"  # Default to Arabidopsis

    if reaction_record:
        reaction_name = reaction_record.get('COMMON-NAME', reaction_id)

        # Get EC number from DBLINKS
        dblinks = reaction_record.get('DBLINKS', [])
        if not isinstance(dblinks, list):
            dblinks = [dblinks] if dblinks else []

        for dblink in dblinks:
            if isinstance(dblink, tuple) and len(dblink) == 2:
                db_name, db_id = dblink
                if db_name == 'EC':
                    ec_number = db_id
                    break

        # Try to get organism info
        species = reaction_record.get('SPECIES', '')
        if species:
            organism = species

    # Create pathway title
    if ec_number:
        title = f"{reaction_name} (EC {ec_number})"
    else:
        title = f"{reaction_name}"

    # Create Xref for the reaction
    xref = None
    if ec_number:
        xref = Xref(identifier=ec_number, dataSource="Enzyme Nomenclature")

    # Create a proper pathway object with all required fields
    pathway = Pathway(
        elementId=builder.id_manager.register_id(f"SINGLE-REACTION-{reaction_id}"),
        title=title,
        organism=organism,
        source="PlantCyc",
        version=None,
        license=None,
        xref=xref,
        description=f"Single reaction view for {reaction_id}: {reaction_name}",
        authors=[Author(name="BioCyc")],
        graphics=Graphics(boardWidth=800.0, boardHeight=600.0),
        dataNodes=[],
        interactions=[],
        graphicalLines=[],
        labels=[],
        shapes=[],
        groups=[],
        annotations=[],
        citations=[],
        evidences=[],
        comments=[Comment(value="This GPML file was automatically generated from BioCyc single reaction data.", source="Automated Conversion")],
        properties=[PathwayProperty(key="UniqueID", value=reaction_id)],
        annotationRefs=[],
        citationRefs=[],
        evidenceRefs=[]
    )

    # Use the same component collection logic as normal pathways
    try:
        pathway_components = builder._collect_pathway_components([reaction_id])

        if not pathway_components['reactions']:
            return None

        # Calculate positions
        positions = builder._calculate_component_positions(pathway_components)

        # Create datanodes
        pathway_datanodes, pathway_groups, compound_node_map, placeholder_stats = builder._create_pathway_datanodes(
            pathway_components, positions
        )

        # Create interactions
        pathway_interactions, new_regulator_groups = builder._create_pathway_interactions(
            pathway_components, compound_node_map, pathway_datanodes
        )

        # Assign to pathway
        pathway.dataNodes = pathway_datanodes
        pathway.groups = pathway_groups + new_regulator_groups
        pathway.interactions = pathway_interactions

        # Collect citations
        element_ids = [reaction_id]

        for compound_id in pathway_components['compounds']:
            element_ids.append(compound_id)

        for protein_id in pathway_components['proteins']:
            element_ids.append(protein_id)

        for gene_id in pathway_components['genes']:
            element_ids.append(gene_id)

        pathway.citations = builder.citation_manager.get_all_citations_for_pathway(element_ids)
        builder._update_pathway_board_size(pathway, pathway_datanodes)

        return pathway

    except Exception as e:
        print(f"    Error building reaction {reaction_id}: {str(e)}")
        return None


def build_single_reactions(builder, unused_reactions, output_dir):
    """Build GPML files for all unused reactions."""
    print(f"\nBuilding {len(unused_reactions)} single reaction GPML files...")

    built_count = 0
    failed_count = 0
    failed_reactions = []

    for i, reaction_id in enumerate(unused_reactions, 1):
        if i % 50 == 0:
            print(f"  Progress: {i}/{len(unused_reactions)} ({i*100//len(unused_reactions)}%)")

        try:
            # Build pathway for single reaction
            pathway = build_single_reaction_pathway(builder, reaction_id)

            if pathway is None:
                failed_count += 1
                failed_reactions.append((reaction_id, "No pathway generated"))
                continue

            # Create safe filename
            safe_reaction_id = re.sub(r'[^a-zA-Z0-9_-]', '_', reaction_id)
            output_filename = f"{safe_reaction_id}.gpml"
            output_filepath = os.path.join(output_dir, output_filename)

            # Export to GPML
            builder.export_pathway_to_gpml(pathway, output_filepath)

            built_count += 1

        except Exception as e:
            failed_count += 1
            failed_reactions.append((reaction_id, str(e)))

    return built_count, failed_count, failed_reactions


def main():
    """Main execution function."""
    # Check command line arguments
    if len(sys.argv) < 3:
        print(__doc__)
        print("\nError: Missing required arguments!")
        print("Usage: python build_pathways.py <data_dir> <output_dir> [--include-reactions] [--layout grid|forceatlas2]")
        print("\nOptions:")
        print("  --include-reactions    : Also build single reaction files for unused reactions")
        print("  --layout grid          : Use grid layout (default)")
        print("  --layout forceatlas2   : Use ForceAtlas2 force-directed layout")
        sys.exit(1)

    data_dir = sys.argv[1]
    output_base_dir = sys.argv[2]
    include_reactions = '--include-reactions' in sys.argv

    # Determine layout type
    layout_type = 'grid'  # default
    for i, arg in enumerate(sys.argv):
        if arg == '--layout' and i + 1 < len(sys.argv):
            layout_type = sys.argv[i + 1].lower()
            if layout_type not in ['grid', 'forceatlas2']:
                print(f"Error: Invalid layout type '{layout_type}'. Use 'grid' or 'forceatlas2'")
                sys.exit(1)

    # Validate data directory
    if not os.path.exists(data_dir):
        print(f"Error: Data directory does not exist: {data_dir}")
        sys.exit(1)

    # Required BioCyc data files
    required_files = [
        'compounds.dat',
        'genes.dat',
        'proteins.dat',
        'reactions.dat',
        'pathways.dat',
        'pubs.dat',
        'regulation.dat'
    ]

    # Check for required files
    missing_files = []
    for filename in required_files:
        filepath = os.path.join(data_dir, filename)
        if not os.path.exists(filepath):
            missing_files.append(filename)

    if missing_files:
        print("Error: Missing required BioCyc data files in data directory:")
        for filename in missing_files:
            print(f"  - {filename}")
        sys.exit(1)

    # Create output directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(output_base_dir, f"biocyc_pathways_{timestamp}")
    individual_pathways_dir = os.path.join(output_dir, "individual_pathways")
    os.makedirs(individual_pathways_dir, exist_ok=True)

    # Create reactions directory if needed
    if include_reactions:
        individual_reactions_dir = os.path.join(output_dir, "individual_reactions")
        os.makedirs(individual_reactions_dir, exist_ok=True)

    print("="*60)
    print("BioCyc to WikiPathways GPML Converter")
    print("="*60)
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Include single reactions: {include_reactions}")
    print(f"Layout type: {layout_type}")
    print("="*60 + "\n")

    # Initialize builder with BioCyc data files
    print("Initializing pathway builder...")
    builder = CompletePathwayBuilderWithGenes(
        compounds_file=os.path.join(data_dir, "compounds.dat"),
        genes_file=os.path.join(data_dir, "genes.dat"),
        proteins_file=os.path.join(data_dir, "proteins.dat"),
        reactions_file=os.path.join(data_dir, "reactions.dat"),
        pathways_file=os.path.join(data_dir, "pathways.dat"),
        pubs_file=os.path.join(data_dir, "pubs.dat"),
        regulation_file=os.path.join(data_dir, "regulation.dat")
    )

    # Find all pathways with reactions
    print("\nFinding all pathways...")
    all_pathways = builder.find_all_pathways()
    print(f"Found {len(all_pathways)} pathways with reactions\n")

    # Build individual pathways
    print("Building individual pathways...")
    built_pathways, failed_pathways = build_individual_pathways(
        builder, all_pathways, individual_pathways_dir, layout_type=layout_type
    )

    # Build single reactions if requested
    built_reactions = 0
    failed_reactions_count = 0
    if include_reactions:
        print("\n" + "="*60)
        print("BUILDING SINGLE REACTIONS")
        print("="*60)

        # Store file paths and reaction processor for later use
        reactions_file = os.path.join(data_dir, "reactions.dat")
        pathways_file = os.path.join(data_dir, "pathways.dat")
        builder.reaction_processor = parsing_utils.read_and_parse(reactions_file)

        # Find unused reactions
        print("\nFinding unused reactions...")
        unused_reactions = find_unused_reactions(builder, reactions_file, pathways_file)
        print(f"  Total reactions: {len(builder.reaction_processor.records)}")
        print(f"  Reactions in pathways: {len(all_pathways)}")
        print(f"  Unused reactions: {len(unused_reactions)}")

        # Build single reactions
        if unused_reactions:
            built_reactions, failed_reactions_count, failed_reactions_list = build_single_reactions(
                builder, unused_reactions, individual_reactions_dir
            )

            # Save build report
            report_file = os.path.join(individual_reactions_dir, "BUILD_REPORT.txt")
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("SINGLE REACTION BUILD REPORT\n")
                f.write("="*70 + "\n\n")
                f.write(f"Build Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write(f"Total unused reactions: {len(unused_reactions)}\n")
                f.write(f"Successfully built: {built_reactions}\n")
                f.write(f"Failed: {failed_reactions_count}\n")
                if len(unused_reactions) > 0:
                    f.write(f"Success rate: {built_reactions*100/len(unused_reactions):.1f}%\n\n")

                if failed_reactions_list:
                    f.write("="*70 + "\n")
                    f.write("FAILED REACTIONS\n")
                    f.write("="*70 + "\n\n")
                    for reaction_id, error in failed_reactions_list:
                        f.write(f"{reaction_id}\n")
                        f.write(f"  Error: {error}\n\n")

    # Print summary to console
    print("\n" + "="*60)
    print("BUILD COMPLETE")
    print("="*60)
    print(f"Pathways built: {len(built_pathways)}/{len(all_pathways)}")
    if include_reactions:
        print(f"Single reactions built: {built_reactions}")
    print(f"Output: {output_dir}")
    print("="*60)

    # Run analysis script automatically
    print("\n" + "="*60)
    print("RUNNING GPML STATISTICS ANALYSIS")
    print("="*60)
    try:
        import subprocess

        print("\nRunning GPML statistics analysis...")
        result = subprocess.run(
            ["python", "analyze_gpml_stats.py", output_dir],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            capture_output=True,
            text=True,
            timeout=300
        )
        if result.returncode == 0:
            print("[+] Analysis completed successfully")
        else:
            print("[-] Analysis failed:")
            if result.stderr:
                print(result.stderr)

    except Exception as e:
        print("[-] Failed to run analysis script: " + str(e))

    if failed_pathways:
        print(f"\nWarning: {len(failed_pathways)} pathway(s) failed to build:")
        for failure in failed_pathways:
            print(f"  - {failure['pathway_id']}: {failure['error']}")

    if include_reactions and failed_reactions_count > 0:
        print(f"\nWarning: {failed_reactions_count} reaction(s) failed to build")
        print(f"See {os.path.join(individual_reactions_dir, 'BUILD_REPORT.txt')} for details")


if __name__ == "__main__":
    main()
