"""
Standard Graphics Styles for PathVisio/WikiPathways

This module provides standard graphics configurations that match
the PathVisio/WikiPathways visual conventions.
"""

from scripts.data_structure.wiki_data_structure import (
    Graphics, HAlign, VAlign, BorderStyle, ShapeType, LineStyle, ConnectorType
)


# Standard DataNode dimensions
STANDARD_WIDTH = 90.0
STANDARD_HEIGHT = 25.0

# Standard font settings
STANDARD_FONT_SIZE = 12.0
STANDARD_FONT_NAME = "Arial"

# Standard colors (hex format without #)
COLOR_BLACK = "000000"
COLOR_WHITE = "FFFFFF"
COLOR_METABOLITE_BLUE = "3955e7"
COLOR_COMPLEX_PURPLE = "9900cc"


def create_standard_datanode_graphics(
    center_x: float,
    center_y: float,
    width: float = STANDARD_WIDTH,
    height: float = STANDARD_HEIGHT,
    color: str = None
) -> Graphics:
    """
    Create standard PathVisio-style DataNode graphics.

    Args:
        center_x: X coordinate of center
        center_y: Y coordinate of center
        width: Width of the node (default: 90.0)
        height: Height of the node (default: 25.0)
        color: Optional fill color for special node types (hex without #)

    Returns:
        Graphics object with standard settings
    """
    return Graphics(
        centerX=center_x,
        centerY=center_y,
        width=width,
        height=height,
        textColor=COLOR_BLACK,
        fontName=STANDARD_FONT_NAME,
        fontWeight=False,
        fontStyle=False,
        fontDecoration=False,
        fontStrikethru=False,
        fontSize=STANDARD_FONT_SIZE,
        hAlign=HAlign.CENTER,
        vAlign=VAlign.MIDDLE,
        borderColor=COLOR_BLACK,
        borderStyle=BorderStyle.SOLID,
        borderWidth=1.0,
        fillColor=color if color else COLOR_WHITE,
        shapeType=ShapeType.RECTANGLE
    )


def create_gene_graphics(center_x: float, center_y: float) -> Graphics:
    """
    Create standard graphics for gene/GeneProduct nodes.

    Args:
        center_x: X coordinate of center
        center_y: Y coordinate of center

    Returns:
        Graphics object with standard gene styling
    """
    return create_standard_datanode_graphics(center_x, center_y)


def create_protein_graphics(center_x: float, center_y: float) -> Graphics:
    """
    Create standard graphics for protein nodes.

    Args:
        center_x: X coordinate of center
        center_y: Y coordinate of center

    Returns:
        Graphics object with standard protein styling
    """
    return create_standard_datanode_graphics(center_x, center_y)


def create_metabolite_graphics(center_x: float, center_y: float) -> Graphics:
    """
    Create standard graphics for metabolite/compound nodes.

    Args:
        center_x: X coordinate of center
        center_y: Y coordinate of center

    Returns:
        Graphics object with standard metabolite styling (blue text, white fill)
    """
    graphics = create_standard_datanode_graphics(center_x, center_y)
    graphics.textColor = COLOR_METABOLITE_BLUE
    graphics.fillColor = COLOR_WHITE
    graphics.borderColor = COLOR_METABOLITE_BLUE
    return graphics


def create_dna_graphics(center_x: float, center_y: float) -> Graphics:
    """
    Create standard graphics for DNA nodes.

    Args:
        center_x: X coordinate of center
        center_y: Y coordinate of center

    Returns:
        Graphics object with standard DNA styling
    """
    return create_standard_datanode_graphics(center_x, center_y)


def create_complex_group_graphics(
    center_x: float,
    center_y: float,
    width: float,
    height: float
) -> Graphics:
    """
    Create standard graphics for protein complex groups.

    Args:
        center_x: X coordinate of center
        center_y: Y coordinate of center
        width: Width of the group
        height: Height of the group

    Returns:
        Graphics object with standard complex group styling
    """
    return Graphics(
        centerX=center_x,
        centerY=center_y,
        width=width,
        height=height,
        textColor="666666",  # Gray text for group labels
        fontName=STANDARD_FONT_NAME,
        fontWeight=False,
        fontStyle=False,
        fontDecoration=False,
        fontStrikethru=False,
        fontSize=STANDARD_FONT_SIZE,
        hAlign=HAlign.CENTER,
        vAlign=VAlign.MIDDLE,
        borderColor=COLOR_COMPLEX_PURPLE,
        borderStyle=BorderStyle.DASHED,
        borderWidth=2.0,
        fillColor="ffffff00",  # Transparent or very light
        shapeType=ShapeType.RECTANGLE,
        zOrder=-1  # Behind other elements
    )


def create_standard_interaction_graphics() -> Graphics:
    """
    Create standard graphics for ALL interactions.

    Returns:
        Graphics object for interactions
    """
    return Graphics(
        lineColor=COLOR_BLACK,
        lineStyle=LineStyle.SOLID,
        lineWidth=1.0,
        connectorType=ConnectorType.STRAIGHT
    )


def create_conversion_graphics() -> Graphics:
    """
    Create standard graphics for conversion/reaction interactions.
    Uses black color with ArrowHead="mim-conversion"
    """
    return create_standard_interaction_graphics()


def create_catalysis_graphics() -> Graphics:
    """
    Create standard graphics for catalysis interactions.
    Uses black color with ArrowHead="mim-catalysis"
    """
    return create_standard_interaction_graphics()


def create_gene_expression_graphics() -> Graphics:
    """
    Create standard graphics for transcription-translation interactions.
    Uses black color with ArrowHead="mim-transcription-translation"
    """
    return create_standard_interaction_graphics()


def create_inhibition_graphics() -> Graphics:
    """
    Create standard graphics for inhibition interactions.
    Uses black color with ArrowHead="mim-inhibition"
    """
    return create_standard_interaction_graphics()


def create_stimulation_graphics() -> Graphics:
    """
    Create standard graphics for stimulation/activation interactions.
    Uses black color with ArrowHead="mim-stimulation"
    """
    return create_standard_interaction_graphics()
