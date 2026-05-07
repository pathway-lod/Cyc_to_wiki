"""
Build Pathways Script

This script handles the execution of pathway building using the CompletePathwayBuilderWithGenes
class. It builds ALL individual pathways, organizing output files
into timestamped directories.
"""

import os
import re
from datetime import datetime
from scripts.build_functions.pathway_builder_core import CompletePathwayBuilderWithGenes


def create_output_directories():
    """
    Create timestamped output directories for organized file storage.

    Returns:
        tuple: (base_dir, individual_dir) paths
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_output_dir = f"biocyc_pathways_{timestamp}"
    individual_pathways_dir = os.path.join(base_output_dir, "individual_pathways")

    os.makedirs(individual_pathways_dir, exist_ok=True)

    return base_output_dir, individual_pathways_dir


def build_individual_pathways(builder, all_pathways, individual_pathways_dir, layout_type='grid', layout_params=None):
    """
    Build all individual pathways from the pathway list.

    Args:
        builder (CompletePathwayBuilderWithGenes): Builder instance
        all_pathways (list): List of pathway info dictionaries
        individual_pathways_dir (str): Directory to save individual pathway files
        layout_type (str): Type of layout to use ('grid' or 'forceatlas2'). Default: 'grid'
        layout_params (dict): Custom parameters for the layout algorithm

    Returns:
        tuple: (built_pathways list, failed_pathways list)
    """
    built_pathways = []
    failed_pathways = []

    for i, pathway_info in enumerate(all_pathways):
        pathway_id = pathway_info['pathway_id']

        try:
            # Build the pathway with specified layout
            pathway = builder.build_complete_pathway_with_genes(pathway_id)

            # Deduplicate elements before exporting
            pathway = builder.deduplicate_pathway_elements(pathway)

            # Check organism attribute
            organism = pathway.organism if pathway.organism else ""

            # Handle missing organism - set to "cellular organisms" (TAX-131567)
            if not organism or organism.strip() == "":
                pathway.organism = "cellular organisms"
                print(f"  Warning: {pathway_id} missing organism, setting to 'cellular organisms'")

            # Create safe filename (just use pathway ID)
            safe_pathway_id = re.sub(r'[^a-zA-Z0-9_-]', '_', pathway_id)
            output_filename = f"{safe_pathway_id}.gpml"
            output_filepath = os.path.join(individual_pathways_dir, output_filename)

            # Export to file
            builder.export_pathway_to_gpml(pathway, output_filepath)

            # Store results
            built_pathways.append({
                'pathway_id': pathway_id,
                'filename': output_filename,
                'filepath': output_filepath
            })

        except Exception as e:
            failed_pathways.append({'pathway_id': pathway_id, 'error': str(e)})

    return built_pathways, failed_pathways