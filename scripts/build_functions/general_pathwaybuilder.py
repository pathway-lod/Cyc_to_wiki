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
            pathway = builder.build_complete_pathway_with_genes(pathway_id, layout_type=layout_type, layout_params=layout_params)

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




def main():
    """
    Main execution function for building all pathways.

    Returns:
        tuple: (builder, built_pathways, failed_pathways)
    """
    # Create output directories
    base_output_dir, individual_pathways_dir = create_output_directories()

    # Initialize builder with BioCyc data files
    builder = CompletePathwayBuilderWithGenes(
        compounds_file="compounds.dat",
        genes_file="genes.dat",
        proteins_file="proteins.dat",
        reactions_file="reactions.dat",
        pathways_file="pathways.dat",
        pubs_file="pubs.dat",
        regulation_file="regulation.dat"
    )

    # Find all pathways
    all_pathways = builder.find_all_pathways()

    # Build individual pathways
    built_pathways, failed_pathways = build_individual_pathways(
        builder, all_pathways, individual_pathways_dir
    )

    # Print summary
    print(f"\nBuilt {len(built_pathways)}/{len(all_pathways)} pathways")
    print(f"Output: {base_output_dir}")

    if failed_pathways:
        print(f"\nWarning: {len(failed_pathways)} pathway(s) failed to build:")
        for failure in failed_pathways:
            print(f"  - {failure['pathway_id']}: {failure['error']}")

    return builder, built_pathways, failed_pathways


if __name__ == "__main__":
    # Run the complete build process
    builder, built_pathways, failed_pathways = main()