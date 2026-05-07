
"""
Layout utilities for positioning pathway elements.
"""
import math

def calculate_component_positions(pathway_components):
    """
    Calculate positions for pathway components in layered layout.
    All layers use grid layout for better organization with no overlap.
    """
    # Constants for spacing
    GENE_SPACING_X = 120
    GENE_SPACING_Y = 60
    PROTEIN_SPACING_X = 150
    PROTEIN_SPACING_Y = 80
    COMPOUND_SPACING_X = 150
    COMPOUND_SPACING_Y = 100
    LAYER_PADDING = 80

    # Layer 1: Genes (top)
    gene_positions = {}
    gene_cols = max(3, math.ceil(math.sqrt(len(pathway_components['genes'])))) if pathway_components['genes'] else 1
    gene_rows = math.ceil(len(pathway_components['genes']) / gene_cols) if pathway_components['genes'] else 0
    gene_start_y = 50
    for i, gene_id in enumerate(pathway_components['genes']):
        row, col = i // gene_cols, i % gene_cols
        gene_positions[gene_id] = (100 + (col * GENE_SPACING_X), gene_start_y + (row * GENE_SPACING_Y))

    gene_layer_end = gene_start_y + (gene_rows * GENE_SPACING_Y) + 25 if gene_rows > 0 else gene_start_y

    # Layer 2: Proteins (middle)
    protein_positions = {}
    protein_cols = max(3, math.ceil(math.sqrt(len(pathway_components['proteins'])))) if pathway_components['proteins'] else 1
    protein_rows = math.ceil(len(pathway_components['proteins']) / protein_cols) if pathway_components['proteins'] else 0
    protein_start_y = gene_layer_end + LAYER_PADDING
    for i, protein_id in enumerate(pathway_components['proteins']):
        row, col = i // protein_cols, i % protein_cols
        protein_positions[protein_id] = (150 + (col * PROTEIN_SPACING_X), protein_start_y + (row * PROTEIN_SPACING_Y))

    protein_layer_end = protein_start_y + (protein_rows * PROTEIN_SPACING_Y) + 25 if protein_rows > 0 else protein_start_y

    # Layer 3: Compounds (bottom)
    compound_positions = {}
    compound_cols = max(3, math.ceil(math.sqrt(len(pathway_components['compounds'])))) if pathway_components['compounds'] else 1
    compound_start_y = protein_layer_end + LAYER_PADDING
    for i, compound_id in enumerate(pathway_components['compounds']):
        row, col = i // compound_cols, i % compound_cols
        compound_positions[compound_id] = (200 + (col * COMPOUND_SPACING_X), compound_start_y + (row * COMPOUND_SPACING_Y))

    return {'genes': gene_positions, 'proteins': protein_positions, 'compounds': compound_positions}

def position_complex_elements(complex_group, monomer_nodes, center_x, center_y):
    """
    Position monomers horizontally within a complex and size the complex group to fit.
    """
    num_monomers = len(monomer_nodes)
    if num_monomers == 0:
        return

    spacing = 90
    total_width = spacing * (num_monomers - 1) + 100

    for i, monomer_node in enumerate(monomer_nodes):
        offset_x = (i - (num_monomers - 1) / 2) * spacing
        monomer_node.graphics.centerX = center_x + offset_x
        monomer_node.graphics.centerY = center_y

    # Update complex group graphics
    group_padding = 20
    complex_group.graphics.centerX = center_x
    complex_group.graphics.centerY = center_y
    complex_group.graphics.width = total_width + group_padding * 2
    complex_group.graphics.height = 60 + group_padding * 2

def update_group_bounds(group, member_nodes):
    """
    Size a group to bound all its member nodes.
    """
    xs = [n.graphics.centerX for n in member_nodes if hasattr(n, 'graphics')]
    ys = [n.graphics.centerY for n in member_nodes if hasattr(n, 'graphics')]
    if not xs or not ys:
        return

    group_padding = 20
    group.graphics.centerX = sum(xs) / len(xs)
    group.graphics.centerY = sum(ys) / len(ys)
    group.graphics.width = (max(xs) - min(xs)) + 100 + group_padding * 2
    group.graphics.height = (max(ys) - min(ys)) + 60 + group_padding * 2
