import uuid
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional

# --- GPML2021 Enums) ---
class DataNodeType(Enum):
    UNDEFINED = "Undefined"
    GENE_PRODUCT = "GeneProduct"
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "Protein"
    COMPLEX = "Complex"
    METABOLITE = "Metabolite"
    PATHWAY = "Pathway"
    DISEASE = "Disease"
    PHENOTYPE = "Phenotype"
    ALIAS = "Alias"
    EVENT = "Event"

class StateType(Enum):
    UNDEFINED = "Undefined"
    PROTEIN_MODIFICATION = "ProteinModification"
    GENETIC_VARIANT = "GeneticVariant"
    EPIGENETIC_MODIFICATION = "EpigeneticModification"

class GroupType(Enum):
    GROUP = "Group"
    TRANSPARENT = "Transparent"
    COMPLEX = "Complex"
    PATHWAY = "Pathway"
    ANALOG = "Analog"
    PARALOG = "Paralog"

class AnnotationType(Enum):
    UNDEFINED = "Undefined"
    ONTOLOGY = "Ontology"
    TAXONOMY = "Taxonomy"

class ArrowHeadType(Enum):
    UNDIRECTED = "Undirected"
    DIRECTED = "Directed"
    CONVERSION = "Conversion"
    INHIBITION = "Inhibition"
    CATALYSIS = "Catalysis"
    STIMULATION = "Stimulation"
    BINDING = "Binding"
    TRANSLOCATION = "Translocation"
    TRANSCRIPTION_TRANSLATION = "TranscriptionTranslation"

class AnchorShapeType(Enum):
    SQUARE = "Square"
    CIRCLE = "Circle"
    NONE = "None"

class HAlign(Enum):
    LEFT = "Left"
    CENTER = "Center"
    RIGHT = "Right"

class VAlign(Enum):
    TOP = "Top"
    MIDDLE = "Middle"
    BOTTOM = "Bottom"

class BorderStyle(Enum):
    SOLID = "Solid"
    DASHED = "Dashed"
    DOUBLE = "Double"

class LineStyle(Enum):
    SOLID = "Solid"
    DASHED = "Dashed"
    DOUBLE = "Double"

class ConnectorType(Enum):
    STRAIGHT = "Straight"
    ELBOW = "Elbow"
    CURVED = "Curved"
    SEGMENTED = "Segmented"

class ShapeType(Enum):
    RECTANGLE = "Rectangle"
    OVAL = "Oval"
    TRIANGLE = "Triangle"
    PENTAGON = "Pentagon"
    HEXAGON = "Hexagon"

# --- Core GPML2021 Data Structures ---

@dataclass
class Xref:
    identifier: str
    dataSource: str

@dataclass
class Url:
    link: str

@dataclass
class Point:
    x: float
    y: float
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    arrowHead: ArrowHeadType = ArrowHeadType.UNDIRECTED
    elementRef: Optional[str] = None
    relX: Optional[float] = None
    relY: Optional[float] = None
@dataclass
@dataclass
class Anchor:
    position: float
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    shapeType: AnchorShapeType = AnchorShapeType.SQUARE

@dataclass
class Author:
    name: str
    username: Optional[str] = None
    order: Optional[int] = None
    xref: Optional[Xref] = None

@dataclass
class Comment:
    value: str
    source: Optional[str] = None

@dataclass
class Property:
    key: str
    value: str

# --- Annotations and Evidence ---
@dataclass
class Annotation:
    value: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    type: AnnotationType = AnnotationType.UNDEFINED
    xref: Optional[Xref] = None
    url: Optional[Url] = None

@dataclass
class Citation:
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    xref: Optional[Xref] = None
    url: Optional[Url] = None

    def __post_init__(self):
        if self.xref is None and self.url is None:
            raise ValueError("Citation must have at least one of 'xref' or 'url'.")

@dataclass
class Evidence:
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    xref: Optional[Xref] = None
    url: Optional[Url] = None
    value: Optional[str] = None

    def __post_init__(self):
        if self.xref is None and self.url is None:
            raise ValueError("Evidence must have at least one of 'xref' or 'url'.")

# --- Reference Structures ---

@dataclass
class AnnotationRef:
    elementRef: str
    citationRefs: List['CitationRef'] = field(default_factory=list)
    evidenceRefs: List['EvidenceRef'] = field(default_factory=list)

@dataclass
class CitationRef:
    elementRef: str
    annotationRefs: List[AnnotationRef] = field(default_factory=list)

@dataclass
class EvidenceRef:
    elementRef: str

# --- Graphics Mixin ---
@dataclass
class Graphics:
    boardWidth: Optional[float] = None
    boardHeight: Optional[float] = None

    centerX: Optional[float] = None
    centerY: Optional[float] = None
    width: Optional[float] = None
    height: Optional[float] = None

    textColor: Optional[str] = None
    fontName: Optional[str] = None
    fontWeight: Optional[bool] = None
    fontStyle: Optional[bool] = None
    fontDecoration: Optional[bool] = None
    fontStrikethru: Optional[bool] = None
    fontSize: Optional[float] = None
    hAlign: Optional[HAlign] = None
    vAlign: Optional[VAlign] = None

    borderColor: Optional[str] = None
    borderStyle: Optional[BorderStyle] = None
    borderWidth: Optional[float] = None
    fillColor: Optional[str] = None
    shapeType: Optional[ShapeType] = None
    zOrder: Optional[int] = None
    rotation: Optional[float] = None

    lineColor: Optional[str] = None
    lineStyle: Optional[LineStyle] = None
    lineWidth: Optional[float] = None
    connectorType: Optional[ConnectorType] = None

# --- CommentGroup Mixin ---

@dataclass
class CommentGroupMixin:
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

# --- Pathway Elements ---

@dataclass
class State:
    textLabel: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    type: StateType = StateType.UNDEFINED
    xref: Optional[Xref] = None
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        req = self.graphics
        if req.relX is None or req.relY is None or req.width is None or req.height is None:
            raise ValueError("State graphics must have 'relX', 'relY', 'width', and 'height'.")

class DataNode(CommentGroupMixin):
    textLabel: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    type: DataNodeType = DataNodeType.UNDEFINED
    groupRef: Optional[str] = None
    aliasRef: Optional[str] = None
    xref: Optional[Xref] = None
    states: List[State] = field(default_factory=list)
    graphics: Graphics = field(default_factory=Graphics)

    def __post_init__(self):
        g = self.graphics
        if g.centerX is None or g.centerY is None or g.width is None or g.height is None:
            raise ValueError("DataNode graphics must have 'centerX', 'centerY', 'width', and 'height'.")

@dataclass
class Interaction:
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    groupRef: Optional[str] = None
    xref: Optional[Xref] = None
    waypoints: List[Point] = field(default_factory=list)
    anchors: List[Anchor] = field(default_factory=list)
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        if len(self.waypoints) < 2:
            raise ValueError("Interaction must have at least 2 waypoints.")
        g = self.graphics
        if g.lineColor is None or g.lineStyle is None or g.lineWidth is None or g.connectorType is None:
            raise ValueError("Interaction graphics must specify 'lineColor', 'lineStyle', 'lineWidth', and 'connectorType'.")

class GraphicalLine:
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    groupRef: Optional[str] = None
    waypoints: List[Point] = field(default_factory=list)
    anchors: List[Anchor] = field(default_factory=list)
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

class Label:
    textLabel: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    groupRef: Optional[str] = None
    href: Optional[str] = None
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        g = self.graphics
        if g.centerX is None or g.centerY is None or g.width is None or g.height is None:
            raise ValueError("Label graphics must have 'centerX', 'centerY', 'width', and 'height'.")

class Shape:
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    textLabel: Optional[str] = None
    groupRef: Optional[str] = None
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        g = self.graphics
        if g.centerX is None or g.centerY is None or g.width is None or g.height is None or g.shapeType is None:
            raise ValueError("Shape graphics must have 'centerX', 'centerY', 'width', 'height', and 'shapeType'.")
@dataclass
class Group:
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    textLabel: Optional[str] = None
    type: GroupType = GroupType.GROUP
    groupRef: Optional[str] = None
    xref: Optional[Xref] = None
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        g = self.graphics
        if g.centerX is None or g.centerY is None or g.width is None or g.height is None:
            raise ValueError("Group graphics must have 'centerX', 'centerY', 'width', and 'height'.")

@dataclass
class Anchor:
    position: float
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    shapeType: AnchorShapeType = AnchorShapeType.SQUARE

@dataclass
class Annotation:
    value: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    type: AnnotationType = AnnotationType.UNDEFINED
    xref: Optional[Xref] = None
    url: Optional[Url] = None

@dataclass
class State:
    textLabel: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    type: StateType = StateType.UNDEFINED
    xref: Optional[Xref] = None
    graphics: Graphics = field(default_factory=Graphics)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        req = self.graphics
        if req.relX is None or req.relY is None or req.width is None or req.height is None:
            raise ValueError("State graphics must have 'relX', 'relY', 'width', and 'height'.")

@dataclass
class DataNode:
    textLabel: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    type: DataNodeType = DataNodeType.UNDEFINED
    groupRef: Optional[str] = None
    aliasRef: Optional[str] = None
    xref: Optional[Xref] = None
    states: List['State'] = field(default_factory=list)
    graphics: 'Graphics' = field(default_factory=Graphics)
    comments: List['Comment'] = field(default_factory=list)
    properties: List['Property'] = field(default_factory=list)
    annotationRefs: List['AnnotationRef'] = field(default_factory=list)
    citationRefs: List['CitationRef'] = field(default_factory=list)
    evidenceRefs: List['EvidenceRef'] = field(default_factory=list)

    def __post_init__(self):
        g = self.graphics
        if g.centerX is None or g.centerY is None or g.width is None or g.height is None:
            raise ValueError("DataNode graphics must have 'centerX', 'centerY', 'width', and 'height'.")

    def __post_init__(self):
        g = self.graphics
        if g.centerX is None or g.centerY is None or g.width is None or g.height is None:
            raise ValueError("DataNode graphics must have 'centerX', 'centerY', 'width', and 'height'.")


@dataclass
class Pathway:
    title: str
    elementId: str = field(default_factory=lambda: str(uuid.uuid4()))
    organism: Optional[str] = None
    source: Optional[str] = None
    version: Optional[str] = None
    license: Optional[str] = None
    xref: Optional[Xref] = None
    description: Optional[str] = None
    authors: List[Author] = field(default_factory=list)
    graphics: Graphics = field(default_factory=Graphics)
    dataNodes: List[DataNode] = field(default_factory=list)
    interactions: List[Interaction] = field(default_factory=list)
    graphicalLines: List[GraphicalLine] = field(default_factory=list)
    labels: List[Label] = field(default_factory=list)
    shapes: List[Shape] = field(default_factory=list)
    groups: List[Group] = field(default_factory=list)
    annotations: List[Annotation] = field(default_factory=list)
    citations: List[Citation] = field(default_factory=list)
    evidences: List[Evidence] = field(default_factory=list)
    comments: List[Comment] = field(default_factory=list)
    properties: List[Property] = field(default_factory=list)
    annotationRefs: List[AnnotationRef] = field(default_factory=list)
    citationRefs: List[CitationRef] = field(default_factory=list)
    evidenceRefs: List[EvidenceRef] = field(default_factory=list)

    def __post_init__(self):
        if self.graphics.boardWidth is None or self.graphics.boardHeight is None:
            raise ValueError("Pathway graphics must have 'boardWidth' and 'boardHeight'.")
