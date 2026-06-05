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
from pathlib import Path

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
from scripts.validate_plantcyc_input import run_checks, fmt_taxon
from datetime import datetime
import csv



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
    # Pathway-level species is always Viridiplantae (plant kingdom).
    # Per-gene/protein species are stored as individual Taxonomy Annotations on DataNodes.
    organism = "Viridiplantae"

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
            from scripts.data_structure.wiki_data_structure import Citation
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

        # Collect DataNode species annotations and add Viridiplantae pathway annotation
        from scripts.data_structure.wiki_data_structure import Annotation, AnnotationType, Xref, AnnotationRef
        pathway_annotations = []
        referenced_anns = set()
        for node in pathway_datanodes + pathway.groups:
            if hasattr(node, 'annotationRefs') and node.annotationRefs:
                for ref in node.annotationRefs:
                    referenced_anns.add(ref.elementRef)
        for ann_id in referenced_anns:
            if ann_id in builder.annotation_index:
                pathway_annotations.append(builder.annotation_index[ann_id])

        viridiplantae_id = "taxonomy_33090"
        if not any(a.elementId == viridiplantae_id for a in pathway_annotations):
            pathway_annotations.append(Annotation(
                elementId=viridiplantae_id,
                value="Viridiplantae",
                type=AnnotationType.TAXONOMY,
                xref=Xref(identifier="33090", dataSource="NCBI Taxonomy")
            ))
        pathway.annotations = pathway_annotations

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


def _write_run_metadata(output_dir, data_dir, db_version, timestamp, include_reactions):
    """Write run.metadata.txt with build provenance."""
    import subprocess

    def _git(cmd):
        try:
            return subprocess.check_output(cmd, cwd=os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))), stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unknown"

    meta_path = os.path.join(output_dir, "run.metadata.txt")
    git_commit = _git(['git','rev-parse','HEAD'])
    git_branch = _git(['git','rev-parse','--abbrev-ref','HEAD'])
    git_dirty  = str(bool(_git(['git','status','--porcelain']))).lower()
    with open(meta_path, "w", encoding="utf-8") as f:
        f.write(f"run_timestamp: {datetime.now().isoformat()}\n")
        f.write(f"repo: Cyc_to_wiki\n")
        f.write(f"git_branch: {git_branch}\n")
        f.write(f"git_commit: {git_commit}\n")
        f.write(f"git_dirty: {git_dirty}\n")
        f.write(f"plantcyc_data_dir: {os.path.abspath(data_dir)}\n")
        f.write(f"plantcyc_version: {db_version or 'unknown'}\n")
        f.write(f"include_reactions: {include_reactions}\n")
        f.write(f"output_dir: {os.path.abspath(output_dir)}\n")
        f.write(f"python: {sys.version.split()[0]}\n")
        f.write(f"command: {' '.join(sys.argv)}\n")
    print(f"  run.metadata.txt written.")


def _generate_species_coverage(data_dir, output_dir, gpml_dir):
    """Run generate_species_coverage_table.py and write both TSV outputs."""
    import subprocess

    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "utils", "generate_species_coverage_table.py")
    if not os.path.exists(script):
        print(f"  WARNING: {script} not found — skipping species coverage tables.")
        return

    cov_orgid_out  = os.path.join(output_dir, "species_coverage_by_plantcyc_orgid.tsv")
    cov_ncbi_out   = os.path.join(output_dir, "species_coverage_by_ncbi.tsv")
    cov_summary_out= os.path.join(output_dir, "species_coverage_summary.tsv")

    cmd = [
        sys.executable, script,
        "--data-dir",       os.path.abspath(data_dir),
        "--gpml-dir",       os.path.abspath(gpml_dir),
        "--output",         os.path.abspath(cov_orgid_out),
        "--output-by-ncbi", os.path.abspath(cov_ncbi_out),
        "--output-summary", os.path.abspath(cov_summary_out),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.stdout:
            for line in result.stdout.strip().splitlines():
                print(f"  {line}")
        if result.returncode != 0 and result.stderr:
            print(f"  WARNING: species coverage script exited {result.returncode}")
            print(f"  {result.stderr[:300]}")
        else:
            print(f"  species_coverage_by_plantcyc_orgid.tsv  (Table S2 — per ORG-ID)")
            print(f"  species_coverage_by_ncbi.tsv            (Table S1 — per NCBI taxon)")
            print(f"  species_coverage_summary.tsv            (summary + absent taxa)")
    except Exception as e:
        print(f"  WARNING: could not generate species coverage tables: {e}")


def _write_validation_report(output_dir, report, data, cross_species_genes, log_lines):
    """Write VALIDATION_REPORT.txt (human-readable) and VALIDATION_SUMMARY.tsv to output_dir."""
    from scripts.validate_plantcyc_input import (
        _load_pubs, _fmt_citation, fmt_taxon, as_list, parse_dat
    )
    import re as _re

    names = data.get("taxon_names", {})
    protein_names = data["protein_names"]
    gene_names    = data["gene_names"]

    # ── Plain-text log ─────────────────────────────────────────────────────
    txt_path = os.path.join(output_dir, "VALIDATION_REPORT.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("PlantCyc input validation — build log\n")
        f.write("="*70 + "\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        for line in log_lines:
            f.write(line + "\n")
        f.write("\n" + "─"*70 + "\n")
        counts = {l: sum(1 for fi in report.findings if fi.level == l)
                  for l in ("ERROR", "WARNING", "INFO")}
        f.write(f"Summary: {counts['ERROR']} error(s), "
                f"{counts['WARNING']} warning(s), {counts['INFO']} info\n")
        if cross_species_genes:
            f.write(f"\nGenes with taxonomy SKIPPED (cross-species products, requires manual curation):\n")
            for g in sorted(cross_species_genes):
                f.write(f"  {g} ({gene_names.get(g, g)})\n")

    # ── TSV supplementary table (S3) ────────────────────────────────────────
    tsv_path = os.path.join(output_dir, "VALIDATION_SUMMARY.tsv")
    HEADER = [
        "check_id", "severity",
        "gene_id", "gene_symbol",
        "protein_id", "protein_name",
        "species_1_taxid", "species_1_name",
        "species_2_taxid", "species_2_name",
        "taxonomy_skipped_in_build",
        "short_note",
        "publication_pmid", "publication_reference",
    ]

    def _split_taxon(tok):
        tok = tok.strip().strip("'\"")
        m = _re.match(r"((?:TAX|ORG)-\S+)\s*\((.+)\)", tok)
        if m:
            return m.group(1), m.group(2)
        return tok, names.get(tok, "")

    # Load citations
    data_dir_path = None
    for fi in report.findings:
        if fi.check in ("gene-cross-species-products",):
            break
    # try to find pubs from already-loaded data
    pubs = {}
    protein_citations = {}

    buf = csv.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(HEADER)

    for finding in report.findings:
        if finding.check not in ("gene-cross-species-products",
                                  "protein-multi-species",
                                  "gene-multi-orgid-products",
                                  "cplx-only-species"):
            continue
        severity = finding.level

        for line in finding.details:
            if finding.check in ("gene-cross-species-products", "gene-multi-orgid-products"):
                m = _re.match(r"(\S+) \(([^)]+)\): (?:products|org_ids)=(\[.*?\]) → (?:species|org_ids)=(\[.*?\])", line)
                if not m:
                    continue
                gid, sym = m.group(1), m.group(2)
                prods = [p.strip().strip("'\"") for p in m.group(3).strip("[]").split(",")]
                sps_raw = m.group(4).strip("[]")
                sp_items = [s.strip() for s in _re.split(r",\s*(?='|TAX|ORG)", sps_raw)]
                tid1, sn1 = _split_taxon(sp_items[0]) if sp_items else ("", "")
                tid2, sn2 = _split_taxon(sp_items[1]) if len(sp_items) > 1 else ("", "")
                prods_str  = "; ".join(prods)
                pnames_str = "; ".join(protein_names.get(p, p) for p in prods)
                skipped = "yes" if gid in cross_species_genes else "no"
                writer.writerow([finding.check, severity, gid, sym,
                                  prods_str, pnames_str, tid1, sn1, tid2, sn2,
                                  skipped, "", "", ""])

            elif finding.check == "protein-multi-species":
                m = _re.match(r"(\S+) \(([^)]+)\): (\[.*)", line)
                if not m:
                    continue
                pid, pname = m.group(1), m.group(2)
                sps_raw = m.group(3)
                sp_items = [s.strip().strip("'\"") for s in
                            _re.split(r",\s*'", sps_raw.strip("[]"))]
                tid1, sn1 = _split_taxon(sp_items[0]) if sp_items else ("", "")
                tid2, sn2 = _split_taxon(sp_items[1]) if len(sp_items) > 1 else ("", "")
                writer.writerow([finding.check, severity, "", "",
                                  pid, pname, tid1, sn1, tid2, sn2,
                                  "no", "", "", ""])

            elif finding.check == "cplx-only-species":
                m = _re.match(r"(\S+) \(([^)]+)\):", line)
                if not m:
                    continue
                taxon_id, taxon_name = m.group(1), m.group(2)
                writer.writerow([finding.check, severity, "", "",
                                  "", "", taxon_id, taxon_name, "", "",
                                  "no",
                                  "Species only annotated via CPLX Group elements",
                                  "", ""])

    with open(tsv_path, "wb") as f:
        f.write(b'\xef\xbb\xbf')  # UTF-8 BOM for Excel
        f.write(buf.getvalue().encode("utf-8"))

    print(f"\n  Validation report: {txt_path}")
    print(f"  Validation TSV:    {tsv_path}")


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
    no_timestamp_subdir = '--no-timestamp-subdir' in sys.argv or '--no_timestamp_subdir' in sys.argv

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

    # Create output directory (optionally without timestamp subdir)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    # Output directory naming:
    #   With --db-version : plantcyc{version}-gpml2021__git{sha7}__{timestamp}
    #   Without           : biocyc_pathways_{timestamp}  (generic fallback)
    # Pass --no-timestamp-subdir (e.g. from run_pipeline.sh) to use output_base_dir directly.
    if no_timestamp_subdir:
        output_dir = output_base_dir
    elif db_version:
        import subprocess as _sp
        try:
            short_sha = _sp.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                stderr=_sp.DEVNULL
            ).decode().strip()
        except Exception:
            short_sha = "unknown"
        output_dir = os.path.join(
            output_base_dir,
            f"plantcyc{db_version}-gpml2021__git{short_sha}__{timestamp}"
        )
    else:
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

    # ── Input validation ──────────────────────────────────────────────────────
    print("="*60)
    print("VALIDATING PLANTCYC INPUT DATA")
    print("="*60)
    validation_report, validation_data = run_checks(Path(data_dir))
    validation_report.print()

    # Collect genes that have products spanning different NCBI species.
    # These genes will have their taxonomy annotation skipped (the species
    # is undefined / pipeline-order-dependent) and flagged in the build log.
    cross_species_genes: set[str] = set()
    for finding in validation_report.findings:
        if finding.check == "gene-cross-species-products":
            import re as _re
            for line in finding.details:
                m = _re.match(r"(\S+) ", line)
                if m:
                    cross_species_genes.add(m.group(1))

    if cross_species_genes:
        print(f"\n  ⚠ {len(cross_species_genes)} gene(s) with cross-species products will have "
              f"their taxonomy annotation SKIPPED during build (manual curation required).")
        print(f"  Affected genes: {sorted(cross_species_genes)}\n")

    # Build log: collect all validation messages for the output report
    build_log_lines = []
    for f in validation_report.findings:
        lvl = f.level
        build_log_lines.append(f"[{lvl}] {f.check}: {f.message}")
        for d in f.details:
            build_log_lines.append(f"    {d}")

    # Write validation report immediately (before building), so it is
    # present even when using --pathway-id / --reaction-id early-exit paths.
    _write_validation_report(output_dir, validation_report, validation_data,
                             cross_species_genes, build_log_lines)
    print("="*60)

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
        organism_mapping=all_mapping,
        cross_species_genes=cross_species_genes,
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

    # ── run.metadata.txt ─────────────────────────────────────────────────────
    _write_run_metadata(output_dir, data_dir, db_version, timestamp, include_reactions)

    # ── Species coverage tables (S1 + S2) ────────────────────────────────────
    print("\n" + "="*60)
    print("GENERATING SPECIES COVERAGE TABLES")
    print("="*60)
    # Pass output_dir (not just individual_pathways_dir) so generate_species_coverage_table.py
    # scans both individual_pathways/*.gpml AND individual_reactions/*.gpml recursively.
    # Species that only appear in reaction files (no full pathway) are otherwise missed.
    _generate_species_coverage(data_dir, output_dir, output_dir)

    # Print summary to console
    print("\n" + "="*60)
    print("BUILD COMPLETE")
    print("="*60)
    print(f"Pathways built: {len(built_pathways)}/{len(all_pathways)}")
    if include_reactions:
        print(f"Single reactions built: {built_reactions}")
    print(f"Output: {output_dir}")
    n_err  = sum(1 for f in validation_report.findings if f.level == "ERROR")
    n_warn = sum(1 for f in validation_report.findings if f.level == "WARNING")
    n_info = sum(1 for f in validation_report.findings if f.level == "INFO")
    print(f"Validation: {n_err} error(s), {n_warn} warning(s), {n_info} info")
    if cross_species_genes:
        print(f"  ⚠ {len(cross_species_genes)} gene(s) had taxonomy skipped (cross-species products)")
    print(f"Output files:")
    print(f"  run.metadata.txt                       GPML_STATISTICS_REPORT.txt")
    print(f"  VALIDATION_REPORT.txt                  VALIDATION_SUMMARY.tsv")
    print(f"  species_coverage_by_ncbi.tsv           species_coverage_summary.tsv")
    print(f"  species_coverage_by_plantcyc_orgid.tsv")
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
