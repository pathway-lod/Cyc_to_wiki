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
    --include-reactions        : Also build single reaction files for unused reactions
    --layout grid              : Use grid layout (default)
    --pathway-id <PATHWAY_ID>  : Build only a specific pathway by ID
    --reaction-id <REACTION_ID>: Build only a specific reaction by ID
    --db-version <VERSION>     : PlantCyc database version (e.g. 17.0.0)

Examples:
    # Build all pathways
    python build_pathways.py ./data ./output

    # Build with DB version
    python build_pathways.py ./data ./output --db-version 17.0.0

    # Build a specific pathway
    python build_pathways.py ./data ./output --pathway-id GLYCOLYSIS

    # Build a specific reaction
    python build_pathways.py ./data ./output --reaction-id RXN-12345

    # Build all with reactions included
    python build_pathways.py ./data ./output --include-reactions

"""

import sys
import io
import os
import re

# Add project root to sys.path to allow imports from scripts.*
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Fix UTF-8 encoding issues on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from scripts.build_functions.general_pathwaybuilder import build_individual_pathways
from scripts.build_functions.pathway_builder_core import CompletePathwayBuilderWithGenes
from scripts.data_structure.wiki_data_structure import Pathway, Graphics, Xref, Author, Comment, Property as PathwayProperty
from scripts.parsing_functions import parsing_utils
from scripts.utils.layout import calculate_component_positions
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
    from scripts.build_functions.build_pathway_data_nodes import convert_to_latin_name

    # Get reaction information
    reaction_record = get_reaction_info(builder, reaction_id)

    # Get reaction name and EC number
    reaction_name = reaction_id
    ec_number = None
    organism = " "

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

        # Try to get organism info and convert to Latin name
        species = reaction_record.get('SPECIES', '')
        if species:
            organism = convert_to_latin_name(str(species))

    # Handle missing organism - set to "cellular organisms"
    if not organism or organism.strip() == "":
        organism = "cellular organisms"

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
        version=builder.version,
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
        positions = calculate_component_positions(pathway_components)

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

        # Add missing citation detection code (same as in build_pathway())
        # Collect all CitationRefs that are used in the pathway
        cited_refs_in_pathway = set()

        for datanode in pathway_datanodes:
            if hasattr(datanode, 'citationRefs') and datanode.citationRefs:
                for ref in datanode.citationRefs:
                    cited_refs_in_pathway.add(ref.elementRef)

        for group in pathway.groups:
            if hasattr(group, 'citationRefs') and group.citationRefs:
                for ref in group.citationRefs:
                    cited_refs_in_pathway.add(ref.elementRef)

        for interaction in pathway.interactions:
            if hasattr(interaction, 'citationRefs') and interaction.citationRefs:
                for ref in interaction.citationRefs:
                    cited_refs_in_pathway.add(ref.elementRef)

        # Add any missing citations that are referenced but not in pathway.citations
        existing_citation_ids = {citation.elementId for citation in pathway.citations if citation.elementId}
        missing_citation_refs = cited_refs_in_pathway - existing_citation_ids

        if missing_citation_refs:
            print(f"  Found {len(missing_citation_refs)} missing citations for reaction {reaction_id}, adding them...")
            from scripts.data_structure.wiki_data_structure import Citation, Xref
            for missing_ref in missing_citation_refs:
                # missing_ref is a sanitized elementId like "citation_PUB_12695547"
                # Check if this citation already exists in citation_objects
                citation = None
                for orig_id, cit_obj in builder.citation_manager.citation_objects.items():
                    if cit_obj.elementId == missing_ref:
                        citation = cit_obj
                        break

                # If not found, try to create it by un-sanitizing the elementId
                if not citation:
                    # Remove "citation_" prefix
                    unsanitized = missing_ref.replace('citation_', '', 1)
                    # Replace underscores with hyphens for BioCyc IDs
                    if unsanitized.startswith('PUB_'):
                        unsanitized = unsanitized.replace('_', '-', 1)  # Only first underscore
                    elif unsanitized.startswith('cit_'):
                        unsanitized = unsanitized.replace('cit_', '', 1)  # Remove cit_ prefix

                    # Try to create citation with un-sanitized ID
                    citation = builder.citation_manager.create_citation_object(unsanitized)

                # Add citation if we got one
                if citation and citation.elementId not in existing_citation_ids:
                    pathway.citations.append(citation)
                    existing_citation_ids.add(citation.elementId)

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

            # Deduplicate elements before exporting
            pathway = builder.deduplicate_pathway_elements(pathway)

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
        print("Usage: python build_pathways.py <data_dir> <output_dir> [options]")
        sys.exit(1)

    data_dir = sys.argv[1]
    output_base_dir = sys.argv[2]
    include_reactions = '--include-reactions' in sys.argv or '--include_reactions' in sys.argv

    # Parse specific pathway or reaction ID
    specific_pathway_id = None
    specific_reaction_id = None
    db_version = None

    for i, arg in enumerate(sys.argv):
        if arg == '--pathway-id' and i + 1 < len(sys.argv):
            specific_pathway_id = sys.argv[i + 1]
        elif arg == '--reaction-id' and i + 1 < len(sys.argv):
            specific_reaction_id = sys.argv[i + 1]
        elif arg == '--db-version' and i + 1 < len(sys.argv):
            db_version = sys.argv[i + 1]

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

    # Attempt to extract version from pathways.dat header if not provided
    if not db_version:
        try:
            pathways_file = os.path.join(data_dir, "pathways.dat")
            if os.path.exists(pathways_file):
                print(f"Attempting to extract version from {pathways_file} header...")
                # Create a temporary reader just for the header
                reader = parsing_utils.FileReader(base_dir=os.path.dirname(data_dir)) # base_dir handling is tricky, giving parent of data_dir
                # Actually FileReader defaults to project root if None.
                # simpler: just read the file directly since we have the full path
                with open(pathways_file, 'r', encoding='latin-1', errors='replace') as f:
                    # Read first few KB to get header
                    head_content = f.read(4096)
                
                header_info = parsing_utils.parse_file_header(head_content)
                if 'Version' in header_info:
                    db_version = header_info['Version']
                    print(f"  Found version: {db_version}")
        except Exception as e:
            print(f"  Warning: Could not extract version from header: {e}")

    # Create output directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
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
    print(f"DB Version: {db_version if db_version else 'Not specified'}")
    print(f"Build Timestamp: {timestamp}")
    print("="*60 + "\n")

    #Build organism mappings from BioCyc classes.dat
    print("="*60)
    print("BUILDING ORGANISM MAPPINGS")
    print("="*60)

    from scripts.utils.build_org_mapping import parse_species_dat

    # Build mapping in memory - {org_id: {'latin_name': str, 'ncbi_id': str}}
    all_mapping = {}

    # Parse classes.dat for both ORG-ID and TAX-ID entries
    classes_dat = os.path.join(data_dir, "classes.dat")
    if os.path.exists(classes_dat):
        print(f"  Reading {classes_dat}...")
        all_mapping.update(parse_species_dat(classes_dat))

    # Also parse species.dat
    species_dat = os.path.join(data_dir, "species.dat")
    if os.path.exists(species_dat):
        print(f"  Reading {species_dat}...")
        all_mapping.update(parse_species_dat(species_dat))

    # Construct final version string
    full_version = f"{db_version}_{timestamp}" if db_version else timestamp

    builder = CompletePathwayBuilderWithGenes(
        compounds_file=os.path.join(data_dir, "compounds.dat"),
        genes_file=os.path.join(data_dir, "genes.dat"),
        proteins_file=os.path.join(data_dir, "proteins.dat"),
        reactions_file=os.path.join(data_dir, "reactions.dat"),
        pathways_file=os.path.join(data_dir, "pathways.dat"),
        pubs_file=os.path.join(data_dir, "pubs.dat"),
        regulation_file=os.path.join(data_dir, "regulation.dat"),
        version=full_version,
        organism_mapping=all_mapping
    )

    # Handle specific pathway build
    if specific_pathway_id:
        print("="*60)
        print(f"BUILDING SPECIFIC PATHWAY: {specific_pathway_id}")
        print("="*60 + "\n")

        # Build only the specific pathway
        # Create the pathway info dict that build_individual_pathways expects
        pathway_list = [{'pathway_id': specific_pathway_id}]
        built_pathways, failed_pathways = build_individual_pathways(
            builder, pathway_list, individual_pathways_dir, layout_type=layout_type
        )

        if built_pathways:
            print(f"\n✓ Successfully built pathway: {specific_pathway_id}")
        else:
            print(f"\n✗ Failed to build pathway: {specific_pathway_id}")
            if failed_pathways:
                print(f"  Error: {failed_pathways[0]['error']}")

        # Skip other build steps
        print("\n" + "="*60)
        print("BUILD COMPLETE")
        print("="*60)
        return

    # Handle specific reaction build
    if specific_reaction_id:
        print("="*60)
        print(f"BUILDING SPECIFIC REACTION: {specific_reaction_id}")
        print("="*60 + "\n")

        # Initialize reaction processor
        reactions_file = os.path.join(data_dir, "reactions.dat")
        builder.reaction_processor = parsing_utils.read_and_parse(reactions_file)

        # Build only the specific reaction
        individual_reactions_dir = os.path.join(output_dir, "individual_reactions")
        os.makedirs(individual_reactions_dir, exist_ok=True)

        pathway = build_single_reaction_pathway(builder, specific_reaction_id)

        if pathway:
            # Deduplicate elements before exporting
            pathway = builder.deduplicate_pathway_elements(pathway)

            # Create safe filename
            safe_reaction_id = re.sub(r'[^a-zA-Z0-9_-]', '_', specific_reaction_id)
            output_filename = f"{safe_reaction_id}.gpml"
            output_filepath = os.path.join(individual_reactions_dir, output_filename)

            # Export to GPML
            builder.export_pathway_to_gpml(pathway, output_filepath)

            print(f"✓ Successfully built reaction: {specific_reaction_id}")
            print(f"  Output: {output_filepath}")
        else:
            print(f"✗ Failed to build reaction: {specific_reaction_id}")
            print(f"  Reaction may not exist or has no data")

        # Skip other build steps
        print("\n" + "="*60)
        print("BUILD COMPLETE")
        print("="*60)
        return

    # Build all pathways (default behavior)
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
            ["python", os.path.join("scripts", "utils", "analyze_gpml_stats.py"), output_dir],
            cwd=project_root,  # Run from project root
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
