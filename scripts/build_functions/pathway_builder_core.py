"""
Core Pathway Builder Module

This module provides the main CompletePathwayBuilderWithGenes class for loading
BioCyc data and building complete pathway networks with genes, proteins, complexes, and reactions.
"""

from scripts.data_structure.wiki_data_structure import (
    Pathway, DataNode, Interaction, Point, Graphics, ArrowHeadType, Anchor, AnchorShapeType,
    LineStyle, ConnectorType, Group, Comment
)
from scripts.build_functions.build_pathway_data_nodes import create_enhanced_pathways_from_file
from scripts.build_functions.build_compounds_data_nodes import create_enhanced_datanodes_from_compounds
from scripts.build_functions.build_protein_data_nodes import create_enhanced_datanodes_from_proteins
from scripts.build_functions.build_gene_data_nodes import create_enhanced_datanodes_from_genes
from scripts.build_functions.build_reaction_interaction import create_all_reaction_interactions
from scripts.build_functions.build_regulation_interaction import create_all_regulation_interactions
from scripts.build_functions.citation_manager import CitationManager
from scripts.object2gmpl.gpml_writer import GPMLWriter
from scripts.parsing_functions import parsing_utils
from scripts.utils import standard_graphics
import math
import re
import copy
import networkx as nx
import os

# Check if forceatlas2 layout is available
try:
    from networkx import forceatlas2_layout
    HAS_FORCEATLAS2 = True
except ImportError:
    HAS_FORCEATLAS2 = False


def sanitize_element_id(element_id):
    """
    Sanitize element IDs for GPML compatibility.

    Args:
        element_id (str): Original element ID

    Returns:
        str: Sanitized element ID safe for GPML
    """
    if not element_id:
        return element_id
    #replace invalid characters with underscores
    sanitized = re.sub(r'[^a-zA-Z0-9_.-]', '_', element_id)
    # check with extra sanitization?
    if sanitized and not (sanitized[0].isalpha() or sanitized[0] == '_'):
        sanitized = f"_{sanitized}"
    return sanitized


class IDManager:
    """
    Manages ID mapping and sanitization for pathway elements.

    Ensures all element IDs are unique and GPML-compatible while maintaining
    a mapping between original and sanitized IDs.
    """

    def __init__(self):
        """Initialize empty ID mapping."""
        self.id_mapping = {}
        self.sanitized_ids = {}  # Track which sanitized IDs are in use: {sanitized_id: original_id}
        self.id_counter = {}  # Counter for handling collisions: {base_id: count}

    def register_id(self, original_id):
        """
        Register and sanitize an ID, ensuring uniqueness.

        Args:
            original_id (str): Original ID from BioCyc data

        Returns:
            str: Sanitized ID safe for GPML (guaranteed unique)
        """
        if original_id not in self.id_mapping:
            sanitized = sanitize_element_id(original_id)

            # Check for collision with existing sanitized IDs
            if sanitized in self.sanitized_ids:
                # Collision detected - append counter to make unique
                base_id = sanitized
                if base_id not in self.id_counter:
                    self.id_counter[base_id] = 1
                else:
                    self.id_counter[base_id] += 1

                # Append counter to make unique
                sanitized = f"{base_id}_{self.id_counter[base_id]}"

            self.id_mapping[original_id] = sanitized
            self.sanitized_ids[sanitized] = original_id

        return self.id_mapping[original_id]

    def get_sanitized_id(self, original_id):
        """
        Get sanitized ID for an original ID.

        Args:
            original_id (str): Original ID to look up

        Returns:
            str: Sanitized ID if found, otherwise original ID
        """
        return self.id_mapping.get(original_id, original_id)

    def get_all_sanitized_ids(self):
        """
        Get all currently registered sanitized IDs.

        Returns:
            set: Set of all sanitized IDs in use
        """
        return set(self.sanitized_ids.keys())


class CompletePathwayBuilderWithGenes:
    """
    Pathway builder for BioCyc data with full gene-protein-reaction networks.

    This class loads all BioCyc data files and creates pathway visualizations
    showing the biological hierarchy: Genes -> Proteins/Complexes -> Reactions -> Compounds.

    """

    def __init__(self, compounds_file, genes_file, proteins_file, reactions_file, pathways_file, pubs_file="pubs.dat", regulation_file="regulation.dat", version=None):
        """
        Initialize builder by loading all BioCyc data files.

        Args:
            compounds_file (str): Path to compounds.dat
            genes_file (str): Path to genes.dat
            proteins_file (str): Path to proteins.dat
            reactions_file (str): Path to reactions.dat
            pathways_file (str): Path to pathways.dat
            pubs_file (str): Path to pubs.dat (default: "pubs.dat")
            regulation_file (str): Path to regulation.dat (default: "regulation.dat")
            version (str): Version string for GPML files (default: None)
        """
        self.id_manager = IDManager()
        self.version = version

        # Initialize CitationManager
        self.citation_manager = CitationManager(pubs_file)
        print(f"Citation Manager initialized with {len(self.citation_manager.publications)} publications")
        
        self.annotation_index = {} # Map elementId -> Annotation object

        # Store file paths
        self.genes_file = genes_file
        self.data_dir = os.path.dirname(genes_file)

        # Load all data types
        self._load_compounds(compounds_file)
        self._load_genes(genes_file)
        self._load_proteins(proteins_file)
        self._load_reactions(reactions_file)
        self._load_regulations(regulation_file)
        self._load_pathways(pathways_file)

        # Build biological mappings
        self._build_mappings()

    def _load_compounds(self, compounds_file):
        """Load and process compound data with citations."""
        self.compound_nodes = create_enhanced_datanodes_from_compounds(compounds_file, self.citation_manager)
        self._register_and_map_nodes(self.compound_nodes)
        self.compound_original_to_node = self._create_original_to_node_map(self.compound_nodes)
        print(f"Loaded {len(self.compound_nodes)} compound_nodes ")

    def _load_genes(self, genes_file):
        """Load and process gene data with citations."""
        self.gene_nodes, self.gene_citations, self.gene_annotations = create_enhanced_datanodes_from_genes(
            genes_file, self.citation_manager
        )
        # Index annotations
        for annotation in self.gene_annotations:
            self.annotation_index[annotation.elementId] = annotation
            
        self._register_and_map_nodes(self.gene_nodes)
        self.gene_original_to_node = self._create_original_to_node_map(self.gene_nodes)

    def _load_proteins(self, proteins_file):
        """Load and process protein data with complex support and enzyme information."""
        # Load proteins with enhanced datanodes and groups for complexes
        result = create_enhanced_datanodes_from_proteins(proteins_file, self.citation_manager)

        if len(result) == 4:
            self.protein_nodes, self.protein_groups, self.protein_citations, self.protein_annotations = result
        elif len(result) == 3:
            # Fallback for old version
            self.protein_nodes, self.protein_groups, self.protein_citations = result
            self.protein_annotations = []
        else:
            # Fallback for even older version
            self.protein_nodes, self.protein_citations = result
            self.protein_groups = []
            self.protein_annotations = []
            
        # Index annotations
        for annotation in self.protein_annotations:
            self.annotation_index[annotation.elementId] = annotation

        print(f"Loaded {len(self.protein_nodes)} protein nodes and {len(self.protein_groups)} complex groups")

        #detect any duplicate elementIds in protein_nodes
        seen_ids = set()
        duplicates = []
        for node in self.protein_nodes:
            if node.elementId in seen_ids:
                duplicates.append(node.elementId)
            seen_ids.add(node.elementId)

        if duplicates:
            print(f"WARNING: Found {len(duplicates)} duplicate protein DataNode IDs: {duplicates}")
            print("Removing duplicates")
            # Remove duplicates by keeping only first occurrence
            unique_nodes = []
            seen_ids = set()
            for node in self.protein_nodes:
                if node.elementId not in seen_ids:
                    unique_nodes.append(node)
                    seen_ids.add(node.elementId)
                else:
                    print(f"  Removing duplicate: {node.elementId}")
            self.protein_nodes = unique_nodes
            print(f"Protein nodes after deduplication: {len(self.protein_nodes)}")

        self._register_and_map_nodes(self.protein_nodes)
        self._register_and_map_nodes(self.protein_groups)

        self.protein_original_to_node = self._create_original_to_node_map(self.protein_nodes)
        self.group_original_to_node = self._create_original_to_node_map(self.protein_groups)

        # Merge complex groups into the protein mapping for enzyme lookups
        for complex_id, complex_group in self.group_original_to_node.items():
            self.protein_original_to_node[complex_id] = complex_group

        # Load raw protein records for enzyme mapping
        self.proteins_processor = parsing_utils.read_and_parse(proteins_file)
        self.protein_records = {r['UNIQUE-ID']: r for r in self.proteins_processor.records if 'UNIQUE-ID' in r}

    def _load_reactions(self, reactions_file):
        """Load and process reaction data with citations."""
        self.reaction_data = create_all_reaction_interactions(reactions_file, self.citation_manager)
        self._register_reaction_ids()
        self.reaction_index = {r['interaction'].elementId: r for r in self.reaction_data}
        print(f"Loaded {len(self.reaction_data )} reactioninteractions")
    def _load_regulations(self, regulation_file):
        """Load and process regulation data with citations."""
        self.regulation_data = create_all_regulation_interactions(regulation_file, self.citation_manager)
        self._register_regulation_ids()
        print(f"Loaded {len(self.regulation_data)} regulation interactions")

    def _load_pathways(self, pathways_file):
        """Load and process pathway data with citations."""
        self.pathways = create_enhanced_pathways_from_file(pathways_file, self.citation_manager, self.version)
        self._register_and_map_nodes(self.pathways)
        self.pathways_processor = parsing_utils.read_and_parse(pathways_file)
        self.pathway_records = {r['UNIQUE-ID']: r for r in self.pathways_processor.records if 'UNIQUE-ID' in r}

    def _build_mappings(self):
        """Build all biological relationship mappings"""
        print("Building gene->protein mappings...")
        self.gene_protein_mapping = self._build_gene_protein_mapping()

        print("Building protein->reaction mappings...")
        self.protein_reaction_mapping, self.reaction_to_enzymes = self._build_complete_enzyme_mapping()

        print("Building regulation mappings...")
        self._build_regulation_mapping()

        self._build_protein_to_genes_mapping()
        print("Mappings complete!")

    def _build_protein_to_genes_mapping(self):
        """
        Build a reverse mapping from proteins to genes for efficient lookup.
        This avoids nested loops in pathway finding.
        """
        self.protein_to_genes = {}

        for gene_id, connections in self.gene_protein_mapping.items():
            for connection in connections:
                protein_id = connection['protein_id']
                if protein_id not in self.protein_to_genes:
                    self.protein_to_genes[protein_id] = set()
                self.protein_to_genes[protein_id].add(gene_id)

        print(f"Built protein->genes mapping with {len(self.protein_to_genes)} entries")

    def _register_and_map_nodes(self, nodes):
        """
        Register and sanitize node IDs for a list of nodes.

        Args:
            nodes (list): List of DataNode or Group objects to process
        """
        for node in nodes:
            original_id = node.elementId
            sanitized_id = self.id_manager.register_id(original_id)
            node.elementId = sanitized_id

    def _create_original_to_node_map(self, nodes):
        """
        Create mapping from original IDs to node objects.

        Args:
            nodes (list): List of DataNode or Group objects

        Returns:
            dict: Mapping from original ID to node object
        """
        original_to_node = {}
        for node in nodes:
            for orig_id, sanit_id in self.id_manager.id_mapping.items():
                if sanit_id == node.elementId:
                    original_to_node[orig_id] = node
                    break
        return original_to_node

    def _register_reaction_ids(self):
        """Register and sanitize all reaction-related IDs."""
        for reaction_dict in self.reaction_data:
            interaction = reaction_dict['interaction']
            original_id = interaction.elementId
            sanitized_id = self.id_manager.register_id(original_id)
            interaction.elementId = sanitized_id

            for waypoint in interaction.waypoints:
                waypoint.elementId = self.id_manager.register_id(waypoint.elementId)
            for anchor in interaction.anchors:
                anchor.elementId = self.id_manager.register_id(anchor.elementId)

    def _register_regulation_ids(self):
        """Register and sanitize all regulation-related IDs."""
        for regulation_dict in self.regulation_data:
            interaction = regulation_dict['interaction']
            original_id = interaction.elementId
            sanitized_id = self.id_manager.register_id(original_id)
            interaction.elementId = sanitized_id

            for waypoint in interaction.waypoints:
                waypoint.elementId = self.id_manager.register_id(waypoint.elementId)
            for anchor in interaction.anchors:
                anchor.elementId = self.id_manager.register_id(anchor.elementId)

    def _build_gene_protein_mapping(self):
        """
        Build mapping from genes to their protein products.

        Returns:
            dict: Maps gene ID to list of protein connection dictionaries
        """
        gene_protein_map = {}

        # Load raw gene records
        genes_processor = parsing_utils.read_and_parse(self.genes_file)

        # Create gene lookup
        gene_by_id = {}
        for gene_record in genes_processor.records:
            gene_id = gene_record.get('UNIQUE-ID')
            if gene_id:
                gene_by_id[gene_id] = gene_record

        # Map genes to their protein products
        genes_with_products = 0
        for gene_id, gene_record in gene_by_id.items():
            products = gene_record.get('PRODUCT', [])
            if not isinstance(products, list):
                products = [products] if products else []

            if products:
                genes_with_products += 1

            gene_protein_map[gene_id] = []

            for product_id in products:
                # Check if we have both the gene and protein nodes
                if gene_id in self.gene_original_to_node and product_id in self.protein_original_to_node:
                    gene_protein_map[gene_id].append({
                        'protein_id': product_id,
                        'gene_node': self.gene_original_to_node[gene_id],
                        'protein_node': self.protein_original_to_node[product_id]
                    })

        return gene_protein_map

    def _build_complete_enzyme_mapping(self):
        """
        Build complete protein->reaction mapping using enzrxns.dat bridge.
        Tracks both proteins and complexes, storing monomer information.
        """
        protein_to_reactions = {}
        reaction_to_enzymes = {}

        # Load enzrxns.dat to get ENZRXN -> RXN mapping
        enzrxns_file = os.path.join(self.data_dir, "enzrxns.dat")
        enzrxns_processor = parsing_utils.read_and_parse(enzrxns_file)

        enzrxn_to_reaction = {}
        for enzrxn_record in enzrxns_processor.records:
            enzrxn_id = enzrxn_record.get('UNIQUE-ID')
            reaction_id = enzrxn_record.get('REACTION')

            if enzrxn_id and reaction_id:
                if isinstance(reaction_id, list):
                    reaction_id = reaction_id[0] if reaction_id else None

                if reaction_id:
                    enzrxn_to_reaction[enzrxn_id] = reaction_id

        print(f"Built enzrxn->reaction mapping with {len(enzrxn_to_reaction)} entries")

        # Get available reaction IDs from loaded reaction data
        available_reaction_ids = set()
        for reaction_dict in self.reaction_data:
            reaction_id = reaction_dict['interaction'].elementId
            # Also check original IDs
            for orig_id, sanit_id in self.id_manager.id_mapping.items():
                if sanit_id == reaction_id:
                    available_reaction_ids.add(orig_id)
                    break

        print(f"Available reactions: {len(available_reaction_ids)}")

        # Build complete mapping chain using raw protein records
        proteins_with_catalyzes = 0
        complexes_with_catalyzes = 0

        for protein_id, protein_record in self.protein_records.items():
            if 'CATALYZES' in protein_record:
                # Check if this is a complex
                types = protein_record.get('TYPES', [])
                if not isinstance(types, list):
                    types = [types] if types else []

                is_complex = any('COMPLEX' in str(t).upper() for t in types)

                if is_complex:
                    complexes_with_catalyzes += 1
                else:
                    proteins_with_catalyzes += 1

                catalyzes = protein_record['CATALYZES']
                if not isinstance(catalyzes, list):
                    catalyzes = [catalyzes] if catalyzes else []

                protein_reactions = []

                for enzyme_id in catalyzes:
                    enzyme_id = str(enzyme_id).strip()

                    if enzyme_id in enzrxn_to_reaction:
                        target_reaction_id = enzrxn_to_reaction[enzyme_id]

                        if target_reaction_id in available_reaction_ids:
                            protein_reactions.append(target_reaction_id)

                            if target_reaction_id not in reaction_to_enzymes:
                                reaction_to_enzymes[target_reaction_id] = []

                            # For complexes, store the group but also track that it's a complex
                            reaction_to_enzymes[target_reaction_id].append({
                                'protein_id': protein_id,
                                'protein_node': self.protein_original_to_node.get(protein_id),
                                'enzyme_id': enzyme_id,
                                'reaction_id': target_reaction_id,
                                'is_complex': is_complex
                            })

                if protein_reactions:
                    protein_to_reactions[protein_id] = protein_reactions

        print(f"Proteins with CATALYZES: {proteins_with_catalyzes}")
        print(f"Complexes with CATALYZES: {complexes_with_catalyzes}")
        print(f"Total reaction->enzyme mappings: {len(reaction_to_enzymes)}")

        # Store monomer nodes for quick lookup
        self.monomer_by_complex = {}
        for node in self.protein_nodes:
            if hasattr(node, 'groupRef') and node.groupRef:
                complex_id = node.groupRef
                if complex_id not in self.monomer_by_complex:
                    self.monomer_by_complex[complex_id] = []
                self.monomer_by_complex[complex_id].append(node)

        print(f"Complex->monomer mappings: {len(self.monomer_by_complex)}")

        return protein_to_reactions, reaction_to_enzymes

    def _build_regulation_mapping(self):
        """
        Build mapping for regulation interactions.

        Maps:
        - regulated_entity (ENZRXN) -> reaction_id
        - regulator -> target (protein, compound, etc.)

        Handles cases where:
        1. REGULATED-ENTITY is an ENZRXN (maps to reaction via enzrxn_to_reaction)
        2. REGULATOR is a compound, protein, or other entity
        """
        # Load enzrxns.dat to get ENZRXN -> RXN mapping (reuse if available)
        if not hasattr(self, 'enzrxn_to_reaction'):
            enzrxns_file = os.path.join(self.data_dir, "enzrxns.dat")
            enzrxns_processor = parsing_utils.read_and_parse(enzrxns_file)

            self.enzrxn_to_reaction = {}
            for enzrxn_record in enzrxns_processor.records:
                enzrxn_id = enzrxn_record.get('UNIQUE-ID')
                reaction_id = enzrxn_record.get('REACTION')

                if enzrxn_id and reaction_id:
                    if isinstance(reaction_id, list):
                        reaction_id = reaction_id[0] if reaction_id else None

                    if reaction_id:
                        self.enzrxn_to_reaction[enzrxn_id] = reaction_id

        # Build regulation mappings
        self.regulation_by_reaction = {}  # reaction_id -> list of regulation data
        self.regulation_by_regulator = {}  # regulator_id -> list of regulation data
        self.regulation_by_regulated_entity = {}  # regulated_entity_id -> list of regulation data

        regulations_mapped_to_reactions = 0
        regulations_with_compound_regulators = 0
        regulations_with_protein_regulators = 0
        regulations_with_missing_entities = 0

        for reg_data in self.regulation_data:
            regulated_entity = reg_data.get('regulated_entity')
            regulator = reg_data.get('regulator')

            # Track regulation by regulated entity
            if regulated_entity:
                if regulated_entity not in self.regulation_by_regulated_entity:
                    self.regulation_by_regulated_entity[regulated_entity] = []
                self.regulation_by_regulated_entity[regulated_entity].append(reg_data)

                # Try to map regulated entity to reaction
                if regulated_entity in self.enzrxn_to_reaction:
                    reaction_id = self.enzrxn_to_reaction[regulated_entity]
                    if reaction_id not in self.regulation_by_reaction:
                        self.regulation_by_reaction[reaction_id] = []
                    self.regulation_by_reaction[reaction_id].append(reg_data)
                    regulations_mapped_to_reactions += 1

            # Track regulation by regulator
            if regulator:
                if regulator not in self.regulation_by_regulator:
                    self.regulation_by_regulator[regulator] = []
                self.regulation_by_regulator[regulator].append(reg_data)

                # Classify regulator type
                if regulator in self.compound_original_to_node:
                    regulations_with_compound_regulators += 1
                elif regulator in self.protein_original_to_node:
                    regulations_with_protein_regulators += 1

            if not regulated_entity or not regulator:
                regulations_with_missing_entities += 1

        print(f"Regulation mappings built:")
        print(f"  Total regulations: {len(self.regulation_data)}")
        print(f"  Mapped to reactions: {regulations_mapped_to_reactions}")
        print(f"  With compound regulators: {regulations_with_compound_regulators}")
        print(f"  With protein regulators: {regulations_with_protein_regulators}")
        print(f"  With missing entities: {regulations_with_missing_entities}")
        print(f"  Unique reactions regulated: {len(self.regulation_by_reaction)}")
        print(f"  Unique regulators: {len(self.regulation_by_regulator)}")

    def _create_gene_protein_interactions(self, gene_id, protein_connections, interaction_counter):
        """
        Create transcription/translation interactions from genes to proteins.
        """
        interactions = []

        for i, connection in enumerate(protein_connections):
            gene_node = connection['gene_node']
            protein_node = connection['protein_node']

            # Handle complexes vs regular proteins
            target_protein_nodes = []
            if isinstance(protein_node, Group):
                if protein_node.elementId in self.monomer_by_complex:
                    monomers = self.monomer_by_complex[protein_node.elementId]
                    if monomers:
                        target_protein_nodes = monomers
                else:
                    continue
            else:
                target_protein_nodes = [protein_node]

            # Check gene node has graphics
            if not hasattr(gene_node, 'graphics') or gene_node.graphics is None:
                continue

            # Create interaction to each target protein node
            for p_idx, target_protein in enumerate(target_protein_nodes):
                if not hasattr(target_protein, 'graphics') or target_protein.graphics is None:
                    continue

                # Use standard graphics for gene expression
                interaction_graphics = standard_graphics.create_gene_expression_graphics()

                gene_protein_interaction = Interaction(
                    elementId=self.id_manager.register_id(f"gene_protein_{interaction_counter}_{i}_{p_idx}"),
                    waypoints=[
                        Point(
                            elementId=self.id_manager.register_id(f"gene_protein_start_{interaction_counter}_{i}_{p_idx}"),
                            x=gene_node.graphics.centerX + gene_node.graphics.width / 2,
                            y=gene_node.graphics.centerY,
                            arrowHead=ArrowHeadType.UNDIRECTED,
                            elementRef=gene_node.elementId,
                            relX=1.0, relY=0.0
                        ),
                        Point(
                            elementId=self.id_manager.register_id(f"gene_protein_end_{interaction_counter}_{i}_{p_idx}"),
                            x=target_protein.graphics.centerX - target_protein.graphics.width / 2,
                            y=target_protein.graphics.centerY,
                            arrowHead=ArrowHeadType.TRANSCRIPTION_TRANSLATION,
                            elementRef=target_protein.elementId,
                            relX=-1.0, relY=0.0
                        )
                    ],
                    anchors=[],
                    graphics=interaction_graphics,
                    comments=[],
                    properties=[],
                    annotationRefs=[],
                    citationRefs=[],
                    evidenceRefs=[]
                )
                interactions.append(gene_protein_interaction)

        return interactions

    def _create_standard_reaction_with_central_anchor(self, reaction_data, compound_node_map, enzyme_nodes,
                                                      original_reaction_id, positioned_node_map, primary_info=None):
        """
        Create reaction visualization with central anchor and enzyme connections.

        Args:
            reaction_data: Reaction data dictionary
            compound_node_map: Map of compound IDs to positioned compound nodes
            enzyme_nodes: List of enzyme information dictionaries
            original_reaction_id: Original reaction ID (before sanitization)
            positioned_node_map: Map of element IDs to positioned DataNode objects in pathway
            primary_info: Optional dict with 'left_primaries', 'right_primaries', and 'direction' fields
                         Example: {'left_primaries': ['CPD-7035'], 'right_primaries': ['PHENYLACETALDEHYDE'], 'direction': 'R2L'}
                         Primary compounds connect directly to the main reaction line, non-primaries go to anchors.

        Implementation:
            - If primary_info is provided, uses the specified primary compounds for the main reaction line
            - The 'direction' field ('L2R' or 'R2L') determines arrow direction:
              * L2R (default): arrow points from left_primaries -> right_primaries
              * R2L: arrow points from right_primaries -> left_primaries (physiologically reversed)
            - For L2R: left_primaries are arrow start (reactants), right_primaries are arrow end (products)
            - For R2L: right_primaries are arrow start, left_primaries are arrow end
            - All other compounds are added as anchored connections (including additional primaries if multiple exist)
            - If multiple primaries exist, only the first one found goes on the main line; others go to anchors
            - If no primary_info, falls back to old behavior (index 0 compounds on main line, rest on anchors)
        """
        interaction = reaction_data['interaction']
        interactions_to_add = []

        # Ensure the main interaction has Graphics (use standard conversion graphics)
        if not hasattr(interaction, 'graphics') or interaction.graphics is None:
            interaction.graphics = standard_graphics.create_conversion_graphics()

        # Find main reactant and product using primary_info if available
        main_reactant = None
        main_product = None

        # Extract primary compound lists and direction
        left_primaries = []
        right_primaries = []
        direction = 'L2R'  # Default to left-to-right
        if primary_info:
            left_primaries = primary_info.get('left_primaries', [])
            right_primaries = primary_info.get('right_primaries', [])
            direction = primary_info.get('direction', 'L2R')

        # Determine which primaries should be arrow start vs arrow end based on direction
        # Direction indicates physiological flow:
        #   L2R (left-to-right): arrow goes from left_primaries -> right_primaries
        #   R2L (right-to-left): arrow goes from right_primaries -> left_primaries (reversed!)
        if direction == 'R2L':
            # Reverse: arrow should point from right primaries to left primaries
            arrow_start_primaries = right_primaries
            arrow_end_primaries = left_primaries
        else:
            # Normal L2R or default: arrow points from left primaries to right primaries
            arrow_start_primaries = left_primaries
            arrow_end_primaries = right_primaries

        # Find main reactant (first compound from arrow start primaries, or fall back to index 0)
        # arrow_start_primaries may be from right_primaries (products) if direction is R2L
        # So we need to search in BOTH reactants and products lists (i think?) check sanity at some point
        if arrow_start_primaries:
            # Search in both reactants and products for the arrow start compound
            for reactant in reaction_data['reactants']:
                if reactant['compound_id'] in arrow_start_primaries:
                    if reactant['compound_id'] in compound_node_map:
                        main_reactant = compound_node_map[reactant['compound_id']]
                        break

            # If not found in reactants, search in products (for R2L case)
            if not main_reactant:
                for product in reaction_data['products']:
                    if product['compound_id'] in arrow_start_primaries:
                        if product['compound_id'] in compound_node_map:
                            main_reactant = compound_node_map[product['compound_id']]
                            break
        else:
            # Fall back to old behavior: use first reactant
            if reaction_data['reactants']:
                reactant_id = reaction_data['reactants'][0]['compound_id']
                if reactant_id in compound_node_map:
                    main_reactant = compound_node_map[reactant_id]

        # Find main product (first compound from arrow end primaries, or fall back to index 0)
        if arrow_end_primaries:
            # Search in both products and reactants for the arrow end compound
            for product in reaction_data['products']:
                if product['compound_id'] in arrow_end_primaries:
                    if product['compound_id'] in compound_node_map:
                        main_product = compound_node_map[product['compound_id']]
                        break

            # If not found in products, search in reactants (for R2L case)
            if not main_product:
                for reactant in reaction_data['reactants']:
                    if reactant['compound_id'] in arrow_end_primaries:
                        if reactant['compound_id'] in compound_node_map:
                            main_product = compound_node_map[reactant['compound_id']]
                            break
        else:
            # Fall back to old behavior: use first product
            if reaction_data['products']:
                product_id = reaction_data['products'][0]['compound_id']
                if product_id in compound_node_map:
                    main_product = compound_node_map[product_id]

        if not (main_reactant and main_product):
            return interactions_to_add

        # Create main reaction line with central anchor
        start_point = Point(
            elementId=self.id_manager.register_id(f"{interaction.elementId}_start"),
            x=main_reactant.graphics.centerX + main_reactant.graphics.width / 2,
            y=main_reactant.graphics.centerY,
            arrowHead=ArrowHeadType.UNDIRECTED,
            elementRef=main_reactant.elementId,
            relX=1.0, relY=0.0
        )

        end_point = Point(
            elementId=self.id_manager.register_id(f"{interaction.elementId}_end"),
            x=main_product.graphics.centerX - main_product.graphics.width / 2,
            y=main_product.graphics.centerY,
            arrowHead=ArrowHeadType.CONVERSION,
            elementRef=main_product.elementId,
            relX=-1.0, relY=0.0
        )

        central_anchor = Anchor(
            elementId=self.id_manager.register_id(f"{interaction.elementId}_anchor"),
            position=0.5,
            shapeType=AnchorShapeType.CIRCLE
        )

        interaction.waypoints = [start_point, end_point]
        interaction.anchors = [central_anchor]
        interactions_to_add.append(interaction)

        # Connect enzymes to central anchor with catalysis arrows
        for i, enzyme_info in enumerate(enzyme_nodes):
            enzyme_entity = enzyme_info['protein_node']
            is_complex = enzyme_info.get('is_complex', False)
            enzrxn_id = enzyme_info.get('enzyme_id', f'enzyme_{i}')  # Use ENZRXN ID or fallback to index

            # Handle complexes vs regular proteins
            target_enzyme_nodes = []
            if isinstance(enzyme_entity, Group) or is_complex:
                complex_id = enzyme_entity.elementId if isinstance(enzyme_entity, Group) else enzyme_info['protein_id']

                # Get template monomers first
                if complex_id in self.monomer_by_complex:
                    template_monomers = self.monomer_by_complex[complex_id]
                    if template_monomers:
                        # Look up the POSITIONED versions of these monomers in the pathway
                        for template_monomer in template_monomers:
                            positioned_monomer = positioned_node_map.get(template_monomer.elementId)
                            if positioned_monomer:
                                target_enzyme_nodes.append(positioned_monomer)

                        if not target_enzyme_nodes:
                            # No positioned monomers found - skip this enzyme
                            continue
                    else:
                        continue
                else:
                    continue
            else:
                target_enzyme_nodes = [enzyme_entity]

            # Create the catalysis interaction for each target enzyme node
            for e_idx, enzyme_node in enumerate(target_enzyme_nodes):
                if not hasattr(enzyme_node, 'graphics') or enzyme_node.graphics is None:
                    continue

                # Use standard catalysis graphics
                enzyme_graphics = standard_graphics.create_catalysis_graphics()

                # Only add monomer index suffix if there are multiple monomers
                if len(target_enzyme_nodes) > 1:
                    interaction_id = f"{interaction.elementId}_{enzrxn_id}_{e_idx}"
                    start_point_id = f"{interaction.elementId}_{enzrxn_id}_start_{e_idx}"
                    end_point_id = f"{interaction.elementId}_{enzrxn_id}_end_{e_idx}"
                else:
                    interaction_id = f"{interaction.elementId}_{enzrxn_id}"
                    start_point_id = f"{interaction.elementId}_{enzrxn_id}_start"
                    end_point_id = f"{interaction.elementId}_{enzrxn_id}_end"

                enzyme_interaction = Interaction(
                    elementId=self.id_manager.register_id(interaction_id),
                    waypoints=[
                        Point(
                            elementId=self.id_manager.register_id(start_point_id),
                            x=enzyme_node.graphics.centerX,
                            y=enzyme_node.graphics.centerY + enzyme_node.graphics.height / 2,
                            arrowHead=ArrowHeadType.UNDIRECTED,
                            elementRef=enzyme_node.elementId,
                            relX=0.0, relY=1.0
                        ),
                        Point(
                            elementId=self.id_manager.register_id(end_point_id),
                            x=0, y=0,
                            arrowHead=ArrowHeadType.CATALYSIS,
                            elementRef=central_anchor.elementId,
                            relX=0.0, relY=-1.0
                        )
                    ],
                    anchors=[],
                    graphics=enzyme_graphics,
                    comments=[],
                    properties=[],
                    annotationRefs=[],
                    citationRefs=[],
                    evidenceRefs=[]
                )
                interactions_to_add.append(enzyme_interaction)

        # Connect additional reactants and products to central anchor
        reaction_xref = interaction.xref if hasattr(interaction, 'xref') else None
        self._add_additional_reactants_products(reaction_data, compound_node_map, central_anchor, interactions_to_add, reaction_xref, primary_info, main_reactant, main_product)

        return interactions_to_add

    def _add_additional_reactants_products(self, reaction_data, compound_node_map, central_anchor, interactions_to_add, reaction_xref=None, primary_info=None, main_reactant=None, main_product=None):
        """
        Add additional reactants and products to central anchor.
        Uses primary_info to determine which compounds should be on anchors vs main line.

        Args:
            reaction_data: Reaction data dictionary
            compound_node_map: Map of compound IDs to compound nodes
            central_anchor: Central anchor point for the reaction
            interactions_to_add: List to append new interactions to
            reaction_xref: Xref from the main reaction to apply to additional interactions
            primary_info: Optional dict with 'left_primaries' and 'right_primaries' lists
            main_reactant: The main reactant DataNode (on the main reaction line)
            main_product: The main product DataNode (on the main reaction line)
        """
        # Determine which compounds are primary (if info available)
        left_primaries = set()
        right_primaries = set()

        if primary_info:
            left_primaries = set(primary_info.get('left_primaries', []))
            right_primaries = set(primary_info.get('right_primaries', []))

        # Additional reactants: add compounds that are NOT on the main line
        # If primary_info exists: skip the FIRST primary (it's on main line), add all others (including other primaries)
        # Otherwise: use old logic (skip index 0)
        main_reactant_id = main_reactant.elementId if main_reactant else None

        for i, reactant in enumerate(reaction_data['reactants']):
            reactant_id = reactant['compound_id']

            # Skip the compound that's on the main line
            if reactant_id in compound_node_map:
                if compound_node_map[reactant_id].elementId == main_reactant_id:
                    continue

            # For old behavior (no primary_info), skip index 0
            if not primary_info and i == 0:
                continue

            if reactant_id in compound_node_map:
                reactant_node = compound_node_map[reactant_id]
                additional_interaction = Interaction(
                    elementId=self.id_manager.register_id(f"{central_anchor.elementId}_reactant_{i}"),
                    xref=reaction_xref,
                    waypoints=[
                        Point(
                            elementId=self.id_manager.register_id(f"{central_anchor.elementId}_reactant_start_{i}"),
                            x=reactant_node.graphics.centerX + reactant_node.graphics.width / 2,
                            y=reactant_node.graphics.centerY,
                            arrowHead=ArrowHeadType.UNDIRECTED,
                            elementRef=reactant_node.elementId,
                            relX=1.0, relY=0.0
                        ),
                        Point(
                            elementId=self.id_manager.register_id(f"{central_anchor.elementId}_reactant_end_{i}"),
                            x=0, y=0,
                            arrowHead=ArrowHeadType.UNDIRECTED,
                            elementRef=central_anchor.elementId,
                            relX=0.0, relY=0.0
                        )
                    ],
                    graphics=standard_graphics.create_conversion_graphics(),
                    anchors=[],
                    comments=[],
                    properties=[],
                    annotationRefs=[],
                    citationRefs=[],
                    evidenceRefs=[]
                )
                interactions_to_add.append(additional_interaction)

        # Additional products: add compounds that are NOT on the main line
        # If primary_info exists: skip the FIRST primary (it's on main line), add all others (including other primaries)
        # Otherwise: use old logic (skip index 0)
        main_product_id = main_product.elementId if main_product else None

        for i, product in enumerate(reaction_data['products']):
            product_id = product['compound_id']

            # Skip the compound that's on the main line
            if product_id in compound_node_map:
                if compound_node_map[product_id].elementId == main_product_id:
                    continue

            # For old behavior (no primary_info), skip index 0
            if not primary_info and i == 0:
                continue

            if product_id in compound_node_map:
                product_node = compound_node_map[product_id]
                additional_interaction = Interaction(
                    elementId=self.id_manager.register_id(f"{central_anchor.elementId}_product_{i}"),
                    xref=reaction_xref,
                    waypoints=[
                        Point(
                            elementId=self.id_manager.register_id(f"{central_anchor.elementId}_product_start_{i}"),
                            x=0, y=0,
                            arrowHead=ArrowHeadType.UNDIRECTED,
                            elementRef=central_anchor.elementId,
                            relX=0.0, relY=0.0
                        ),
                        Point(
                            elementId=self.id_manager.register_id(f"{central_anchor.elementId}_product_end_{i}"),
                            x=product_node.graphics.centerX - product_node.graphics.width / 2,
                            y=product_node.graphics.centerY,
                            arrowHead=ArrowHeadType.CONVERSION,
                            elementRef=product_node.elementId,
                            relX=-1.0, relY=0.0
                        )
                    ],
                    graphics=standard_graphics.create_conversion_graphics(),
                    anchors=[],
                    comments=[],
                    properties=[],
                    annotationRefs=[],
                    citationRefs=[],
                    evidenceRefs=[]
                )
                interactions_to_add.append(additional_interaction)

    def _parse_primaries(self, primaries_field):
        """
        Parse PRIMARIES field from pathway record.

        Example format: ("RXN-11659" ("CPD-110") ("CPD-12629"))
        Means: RXN-11659 has CPD-110 as primary reactant and CPD-12629 as primary product

        Args:
            primaries_field: List of PRIMARIES entries from pathway record

        Returns:
            dict: Mapping {reaction_id: {'reactants': [...], 'products': [...]}}
        """
        primary_map = {}

        if not primaries_field:
            return primary_map

        # Ensure it's a list
        if not isinstance(primaries_field, list):
            primaries_field = [primaries_field]

        for entry in primaries_field:
            entry_str = str(entry)

            # Parse pattern: ("RXN-ID" ("COMPOUND1" "COMPOUND2") ("COMPOUND3"))
            # First part is reaction ID, second part is reactants, third part is products
            match = re.match(r'\("([^"]+)"\s+\(([^)]+)\)\s+\(([^)]+)\)\)', entry_str)
            if match:
                reaction_id = match.group(1)
                reactants_str = match.group(2)
                products_str = match.group(3)

                # Parse compound IDs (remove quotes and split)
                reactants = [c.strip().strip('"') for c in reactants_str.split() if c.strip().strip('"')]
                products = [c.strip().strip('"') for c in products_str.split() if c.strip().strip('"')]

                primary_map[reaction_id] = {
                    'reactants': reactants,
                    'products': products
                }

        return primary_map

    def _parse_reaction_layout(self, layout_field):
        """
        Parse REACTION-LAYOUT field from pathway record.

        Example format: (RXN-11659 (:LEFT-PRIMARIES CPD-110) (:DIRECTION :L2R) (:RIGHT-PRIMARIES CPD-12629))

        Args:
            layout_field: List of REACTION-LAYOUT entries from pathway record

        Returns:
            dict: Mapping {reaction_id: {'left_primaries': [...], 'right_primaries': [...], 'direction': '...'}}
        """
        layout_map = {}

        if not layout_field:
            return layout_map

        # Ensure it's a list
        if not isinstance(layout_field, list):
            layout_field = [layout_field]

        for entry in layout_field:
            entry_str = str(entry)

            # Parse pattern: (RXN-ID (:LEFT-PRIMARIES CPD1 CPD2) (:DIRECTION :L2R) (:RIGHT-PRIMARIES CPD3))
            reaction_match = re.match(r'\(([^\s]+)', entry_str)
            if not reaction_match:
                continue

            reaction_id = reaction_match.group(1)

            # Extract LEFT-PRIMARIES
            left_match = re.search(r':LEFT-PRIMARIES\s+([^)]+)', entry_str)
            left_primaries = []
            if left_match:
                left_str = left_match.group(1)
                left_primaries = [c.strip() for c in left_str.split() if c.strip() and not c.startswith(':')]

            # Extract RIGHT-PRIMARIES
            right_match = re.search(r':RIGHT-PRIMARIES\s+([^)]+)', entry_str)
            right_primaries = []
            if right_match:
                right_str = right_match.group(1)
                right_primaries = [c.strip() for c in right_str.split() if c.strip() and not c.startswith(':')]

            # Extract DIRECTION
            direction_match = re.search(r':DIRECTION\s+:([^\s)]+)', entry_str)
            direction = direction_match.group(1) if direction_match else 'L2R'

            # Store the layout information
            layout_map[reaction_id] = {
                'left_primaries': left_primaries,
                'right_primaries': right_primaries,
                'direction': direction
            }

        return layout_map

    def _collect_pathway_components(self, reaction_list, pathway_record=None):
        """
        Collect all components (genes, proteins, compounds, reactions) for a pathway.

        Args:
            reaction_list: List of reaction IDs in the pathway
            pathway_record: Optional pathway record containing PRIMARIES and REACTION-LAYOUT
        """
        pathway_compounds = set()
        pathway_reactions = []
        pathway_proteins = set()
        pathway_genes = set()

        # Parse primary compound information from ALL pathway records that might contain these reactions
        # This is important for sub-pathways where REACTION-LAYOUT is defined in the sub-pathway, not the parent
        primary_info = {}

        # First, collect from the provided pathway_record
        if pathway_record:
            primaries_field = pathway_record.get('PRIMARIES', [])
            primaries_map = self._parse_primaries(primaries_field)

            layout_field = pathway_record.get('REACTION-LAYOUT', [])
            layout_map = self._parse_reaction_layout(layout_field)

            # Merge both sources
            for rxn_id in set(list(primaries_map.keys()) + list(layout_map.keys())):
                primary_info[rxn_id] = {
                    'left_primaries': [],
                    'right_primaries': []
                }

                if rxn_id in layout_map:
                    primary_info[rxn_id]['left_primaries'] = layout_map[rxn_id]['left_primaries']
                    primary_info[rxn_id]['right_primaries'] = layout_map[rxn_id]['right_primaries']
                    primary_info[rxn_id]['direction'] = layout_map[rxn_id].get('direction', 'L2R')
                elif rxn_id in primaries_map:
                    primary_info[rxn_id]['left_primaries'] = primaries_map[rxn_id]['reactants']
                    primary_info[rxn_id]['right_primaries'] = primaries_map[rxn_id]['products']

        # This handles the case where a reaction is in a sub-pathway that has its own REACTION-LAYOUT
        # REACTION-LAYOUT is more accurate than PRIMARIES (probably), so we search for it even if we already have PRIMARIES
        for reaction_id in reaction_list:
            # Check if we already have a REACTION-LAYOUT (with direction field) for this reaction
            # If we only have PRIMARIES, we should try to find a better REACTION-LAYOUT
            has_layout = reaction_id in primary_info and 'direction' in primary_info[reaction_id]

            if not has_layout:
                # Search all pathway records for this reaction's layout
                for path_id, path_record in self.pathway_records.items():
                    # Check if this pathway has REACTION-LAYOUT for our reaction
                    layout_field = path_record.get('REACTION-LAYOUT', [])
                    if layout_field:
                        layout_map = self._parse_reaction_layout(layout_field)
                        if reaction_id in layout_map:
                            # REACTION-LAYOUT found - this takes precedence over PRIMARIES
                            primary_info[reaction_id] = {
                                'left_primaries': layout_map[reaction_id]['left_primaries'],
                                'right_primaries': layout_map[reaction_id]['right_primaries'],
                                'direction': layout_map[reaction_id].get('direction', 'L2R')
                            }
                            break

                # If still not found, check PRIMARIES as fallback
                if reaction_id not in primary_info:
                    for path_id, path_record in self.pathway_records.items():
                        primaries_field = path_record.get('PRIMARIES', [])
                        if primaries_field:
                            primaries_map = self._parse_primaries(primaries_field)
                            if reaction_id in primaries_map:
                                primary_info[reaction_id] = {
                                    'left_primaries': primaries_map[reaction_id]['reactants'],
                                    'right_primaries': primaries_map[reaction_id]['products']
                                }
                                break  # Found PRIMARIES, stop searching

        # From reactions -> proteins -> genes
        for original_reaction_id in reaction_list:
            sanitized_reaction_id = self.id_manager.get_sanitized_id(original_reaction_id)

            if sanitized_reaction_id in self.reaction_index:
                reaction_data = self.reaction_index[sanitized_reaction_id]
                pathway_reactions.append((reaction_data, original_reaction_id))

                # Collect compounds
                for reactant in reaction_data['reactants']:
                    pathway_compounds.add(reactant['compound_id'])
                for product in reaction_data['products']:
                    pathway_compounds.add(product['compound_id'])

                # Collect proteins that catalyze this reaction
                if original_reaction_id in self.reaction_to_enzymes:
                    for enzyme_info in self.reaction_to_enzymes[original_reaction_id]:
                        pathway_proteins.add(enzyme_info['protein_id'])

        # From proteins -> genes (using optimized mapping)
        for protein_id in pathway_proteins:
            if protein_id in self.protein_to_genes:
                pathway_genes.update(self.protein_to_genes[protein_id])

        return {
            'genes': pathway_genes,
            'proteins': pathway_proteins,
            'compounds': pathway_compounds,
            'reactions': pathway_reactions,
            'primary_info': primary_info
        }

    def _calculate_component_positions(self, pathway_components):
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
        LAYER_PADDING = 80  # Padding between layers

        # Layer 1: Genes (top) - Grid layout
        gene_positions = {}
        gene_cols = max(3, math.ceil(math.sqrt(len(pathway_components['genes'])))) if pathway_components['genes'] else 1
        gene_rows = math.ceil(len(pathway_components['genes']) / gene_cols) if pathway_components['genes'] else 0

        gene_start_y = 50
        for i, gene_id in enumerate(pathway_components['genes']):
            row = i // gene_cols
            col = i % gene_cols
            x = 100 + (col * GENE_SPACING_X)
            y = gene_start_y + (row * GENE_SPACING_Y)
            gene_positions[gene_id] = (x, y)

        # Calculate where gene layer ends (last row y + node height + padding)
        gene_layer_end = gene_start_y + (gene_rows * GENE_SPACING_Y) + 25 if gene_rows > 0 else gene_start_y

        # Layer 2: Proteins (middle) - Grid layout, starts after gene layer
        protein_positions = {}
        protein_cols = max(3, math.ceil(math.sqrt(len(pathway_components['proteins'])))) if pathway_components['proteins'] else 1
        protein_rows = math.ceil(len(pathway_components['proteins']) / protein_cols) if pathway_components['proteins'] else 0

        protein_start_y = gene_layer_end + LAYER_PADDING
        for i, protein_id in enumerate(pathway_components['proteins']):
            row = i // protein_cols
            col = i % protein_cols
            x = 150 + (col * PROTEIN_SPACING_X)
            y = protein_start_y + (row * PROTEIN_SPACING_Y)
            protein_positions[protein_id] = (x, y)

        # Calculate where protein layer ends
        protein_layer_end = protein_start_y + (protein_rows * PROTEIN_SPACING_Y) + 25 if protein_rows > 0 else protein_start_y

        # Layer 3: Compounds (bottom) - Grid layout, starts after protein layer
        compound_positions = {}
        compound_cols = max(3, math.ceil(math.sqrt(len(pathway_components['compounds'])))) if pathway_components['compounds'] else 1

        compound_start_y = protein_layer_end + LAYER_PADDING
        for i, compound_id in enumerate(pathway_components['compounds']):
            row = i // compound_cols
            col = i % compound_cols
            x = 200 + (col * COMPOUND_SPACING_X)
            y = compound_start_y + (row * COMPOUND_SPACING_Y)
            compound_positions[compound_id] = (x, y)

        return {
            'genes': gene_positions,
            'proteins': protein_positions,
            'compounds': compound_positions
        }

    def _calculate_component_positions_forceatlas2(self, pathway_components, datanodes, layout_params="ForceAtlas2"):
        """
        Calculate positions for pathway components using ForceAtlas2 layout algorithm.

        This method uses force-directed graph layout to automatically position nodes
        based on their connections, resulting in a more organic visualization.

        Args:
            pathway_components (dict): Dict with 'genes', 'proteins', 'compounds', 'reactions'
            datanodes (list): List of DataNode objects for the pathway
            layout_params (dict, optional): Custom ForceAtlas2 parameters

        Returns:
            dict: Positions dict with 'genes', 'proteins', 'compounds' keys mapping IDs to (x, y)
        """
        if not HAS_FORCEATLAS2:
            print("Warning: ForceAtlas2 not available. Falling back to grid layout.")
            return self._calculate_component_positions(pathway_components)

        # Build NetworkX graph from datanodes and their connections
        graph = nx.Graph()

        # Create a mapping of node IDs to track which component type they belong to
        id_to_component = {}

        # Add all nodes to the graph
        for gene_id in pathway_components['genes']:
            graph.add_node(gene_id)
            id_to_component[gene_id] = 'gene'

        for protein_id in pathway_components['proteins']:
            graph.add_node(protein_id)
            id_to_component[protein_id] = 'protein'

        for compound_id in pathway_components['compounds']:
            graph.add_node(compound_id)
            id_to_component[compound_id] = 'compound'

        # Add edges based on reactions and interactions
        # Get all reactions
        for reaction_data, original_reaction_id in pathway_components['reactions']:
            # Extract compounds from the reaction
            left_compounds = reaction_data.get('left_primaries', [])
            right_compounds = reaction_data.get('right_primaries', [])

            # Connect all left compounds to all right compounds
            for left_cpd in left_compounds:
                for right_cpd in right_compounds:
                    if left_cpd in id_to_component and right_cpd in id_to_component:
                        graph.add_edge(left_cpd, right_cpd)

        # If graph is empty or has very few nodes, fall back to grid
        if graph.number_of_nodes() < 2:
            return self._calculate_component_positions(pathway_components)

        # Apply ForceAtlas2 layout
        if layout_params is None:
            layout_params = {
                'max_iter': 500,
                'scaling_ratio': 10.0,
                'gravity': 2.0,
                'jitter_tolerance': 1.0,
                'seed': 42,
            }

        try:
            pos = forceatlas2_layout(graph, **layout_params)
        except Exception as e:
            print(f"Error applying ForceAtlas2 layout: {e}. Falling back to grid layout.")
            return self._calculate_component_positions(pathway_components)

        # Scale positions to reasonable GPML coordinates (0-2000 range)
        if pos:
            x_coords = [p[0] for p in pos.values()]
            y_coords = [p[1] for p in pos.values()]

            x_min, x_max = min(x_coords), max(x_coords)
            y_min, y_max = min(y_coords), max(y_coords)

            x_range = x_max - x_min if x_max != x_min else 1
            y_range = y_max - y_min if y_max != y_min else 1

            # Scale to GPML coordinates
            scale_factor = 500
            scaled_pos = {}
            for node_id, (x, y) in pos.items():
                x_norm = (x - x_min) / x_range
                y_norm = (y - y_min) / y_range

                x_scaled = 100 + (x_norm * scale_factor * 3)
                y_scaled = 50 + (y_norm * scale_factor * 3)

                scaled_pos[node_id] = (x_scaled, y_scaled)

            pos = scaled_pos

        # Organize positions by component type
        gene_positions = {gid: pos.get(gid, (100, 50)) for gid in pathway_components['genes']}
        protein_positions = {pid: pos.get(pid, (150, 150)) for pid in pathway_components['proteins']}
        compound_positions = {cid: pos.get(cid, (200, 250)) for cid in pathway_components['compounds']}

        return {
            'genes': gene_positions,
            'proteins': protein_positions,
            'compounds': compound_positions
        }

    def _create_pathway_datanodes(self, pathway_components, positions):
        """
        Create DataNode objects for pathway components, with Groups as visual overlays.

        Returns:
            tuple: (list of DataNodes, list of Groups, compound node map, placeholder stats)
        """
        import copy
        from scripts.data_structure.wiki_data_structure import HAlign, VAlign, BorderStyle, ShapeType

        pathway_datanodes = []
        pathway_groups = []
        compound_node_map = {}

        # Track placeholders
        placeholder_stats = {
            'compounds': []
        }

        # Add genes
        for gene_id in pathway_components['genes']:
            if gene_id in self.gene_original_to_node:
                gene_node = copy.deepcopy(self.gene_original_to_node[gene_id])

                # Use standard gene graphics if missing
                if not hasattr(gene_node, 'graphics') or gene_node.graphics is None:
                    gene_node.graphics = standard_graphics.create_gene_graphics(0.0, 0.0)

                # Update position
                if gene_id in positions['genes']:
                    x, y = positions['genes'][gene_id]
                    gene_node.graphics.centerX = x
                    gene_node.graphics.centerY = y

                pathway_datanodes.append(gene_node)

        # Process proteins and complexes
        processed_proteins = set()
        global_monomer_ids = set()  # Track ALL monomers added across ALL complexes (prevents duplicates)

        for protein_id in pathway_components['proteins']:
            # Check if this is a complex (either in group mapping or if protein entity is a Group)
            protein_entity = self.protein_original_to_node.get(protein_id)
            is_complex_group = protein_id in self.group_original_to_node or isinstance(protein_entity, Group)

            if is_complex_group:
                # It's a complex - create group and monomers
                if protein_id in self.group_original_to_node:
                    complex_group = copy.deepcopy(self.group_original_to_node[protein_id])
                else:
                    # protein_entity is a Group but not in group mapping
                    complex_group = copy.deepcopy(protein_entity)

                # Get position for this complex
                if protein_id in positions['proteins']:
                    complex_x, complex_y = positions['proteins'][protein_id]
                else:
                    complex_x, complex_y = 150, 150  # Default position

                # Find and position monomer components
                monomer_nodes_in_complex = []
                complex_monomer_ids = set()  # Track monomers for THIS complex

                # Look for monomers with groupRef pointing to this complex
                for node in self.protein_nodes:
                    if hasattr(node, 'groupRef'):
                        # Check both original and sanitized IDs
                        if node.groupRef == complex_group.elementId or \
                           (protein_id in self.id_manager.id_mapping and
                            node.groupRef == self.id_manager.id_mapping[protein_id]):

                            # Skip if we've already added this monomer globally
                            if node.elementId in global_monomer_ids:
                                continue

                            # Skip if we've already added this monomer to this complex
                            if node.elementId in complex_monomer_ids:
                                continue

                            monomer_node = copy.deepcopy(node)

                            # Use standard protein graphics if missing
                            if not hasattr(monomer_node, 'graphics') or monomer_node.graphics is None:
                                monomer_node.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)

                            monomer_nodes_in_complex.append(monomer_node)
                            complex_monomer_ids.add(node.elementId)
                            global_monomer_ids.add(node.elementId)  # Track globally
                            processed_proteins.add(node.elementId)

                # Position monomers within the complex area
                num_monomers = len(monomer_nodes_in_complex)
                if num_monomers > 0:
                    # Arrange monomers horizontally within the complex
                    spacing = 90  # Space between monomer centers
                    total_width = spacing * (num_monomers - 1) + 100

                    for i, monomer_node in enumerate(monomer_nodes_in_complex):
                        # Calculate monomer position
                        offset_x = (i - (num_monomers - 1) / 2) * spacing
                        monomer_node.graphics.centerX = complex_x + offset_x
                        monomer_node.graphics.centerY = complex_y

                        pathway_datanodes.append(monomer_node)

                    # Set group graphics to overlay all monomers
                    group_padding = 20
                    group_width = total_width + group_padding * 2
                    group_height = 60 + group_padding * 2

                    complex_group.graphics = Graphics(
                        centerX=complex_x,
                        centerY=complex_y,
                        width=group_width,
                        height=group_height,
                        textColor='666666',  # Gray text
                        fontName='Arial',
                        fontWeight=False,
                        fontStyle=False,
                        fontDecoration=False,
                        fontStrikethru=False,
                        fontSize=10.0,
                        hAlign=HAlign.CENTER,
                        vAlign=VAlign.MIDDLE,
                        borderColor='9900CC',
                        borderStyle=BorderStyle.DASHED,
                        borderWidth=2.0,
                        fillColor='F5E6FF',
                        shapeType=ShapeType.RECTANGLE,
                        zOrder=-1
                    )

                    pathway_groups.append(complex_group)

            elif protein_id in self.protein_original_to_node:
                # Regular protein (not a complex)
                if protein_id not in processed_proteins:
                    protein_entity = self.protein_original_to_node[protein_id]

                    # Skip if this is actually a Group object (safety check)
                    if isinstance(protein_entity, Group):
                        continue

                    protein_node = copy.deepcopy(protein_entity)

                    # Use standard protein graphics if missing
                    if not hasattr(protein_node, 'graphics') or protein_node.graphics is None:
                        protein_node.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)

                    if protein_id in positions['proteins']:
                        x, y = positions['proteins'][protein_id]
                        protein_node.graphics.centerX = x
                        protein_node.graphics.centerY = y

                    pathway_datanodes.append(protein_node)
                    processed_proteins.add(protein_id)

        # Add compounds (including placeholders for missing ones)
        for compound_id in pathway_components['compounds']:
            if compound_id in self.compound_original_to_node:
                compound_node = copy.deepcopy(self.compound_original_to_node[compound_id])

                # Use standard metabolite graphics if missing
                if not hasattr(compound_node, 'graphics') or compound_node.graphics is None:
                    compound_node.graphics = standard_graphics.create_metabolite_graphics(0.0, 0.0)

                # Update position
                if compound_id in positions['compounds']:
                    x, y = positions['compounds'][compound_id]
                    compound_node.graphics.centerX = x
                    compound_node.graphics.centerY = y

                pathway_datanodes.append(compound_node)
                compound_node_map[compound_id] = compound_node
            else:
                # Create placeholder node for missing compound
                from scripts.data_structure.wiki_data_structure import DataNode, Xref, ShapeType, DataNodeType, Property
                from scripts.utils.HTML_cleaner import clean_text_label

                sanitized_id = self.id_manager.register_id(compound_id)
                placeholder_node = DataNode(
                    elementId=sanitized_id,
                    textLabel=clean_text_label(compound_id),
                    type=DataNodeType.METABOLITE,
                    xref=Xref(identifier=compound_id, dataSource=''),
                    graphics=standard_graphics.create_metabolite_graphics(0.0, 0.0),
                    properties=[Property(key='IsPlaceholder', value='true')]
                )

                # Set visual indication for placeholder (lighter color)
                placeholder_node.graphics.fillColor = 'F0F0F0'  # Light gray
                placeholder_node.graphics.borderColor = 'AAAAAA'  # Gray border

                # Update position
                if compound_id in positions['compounds']:
                    x, y = positions['compounds'][compound_id]
                    placeholder_node.graphics.centerX = x
                    placeholder_node.graphics.centerY = y

                pathway_datanodes.append(placeholder_node)
                compound_node_map[compound_id] = placeholder_node

                # Track placeholder
                placeholder_stats['compounds'].append(compound_id)

        # Check for orphaned groupRefs (monomers that reference groups not in pathway_groups)
        group_ids_in_pathway = {group.elementId for group in pathway_groups}
        datanode_ids_in_pathway = {node.elementId for node in pathway_datanodes}
        missing_group_refs = set()

        # First pass: identify missing groups
        for datanode in pathway_datanodes:
            if hasattr(datanode, 'groupRef') and datanode.groupRef:
                if datanode.groupRef not in group_ids_in_pathway:
                    missing_group_refs.add(datanode.groupRef)

        # Second pass: fix nested complexes (complexes that are both DataNodes and have monomers)
        # For these, remove the groupRef from their monomers to avoid orphaned references
        nested_complexes = missing_group_refs & datanode_ids_in_pathway

        if nested_complexes:
            # Remove groupRef from monomers of nested complexes
            for datanode in pathway_datanodes:
                if hasattr(datanode, 'groupRef') and datanode.groupRef in nested_complexes:
                    datanode.groupRef = None

        # Only create placeholder groups for truly missing ones (not nested complexes)
        missing_group_refs = missing_group_refs - nested_complexes

        # Create placeholder groups for missing groupRefs
        for missing_group_id in missing_group_refs:
            # Try to find the group in our group mappings
            group_template = None

            # Check in group_original_to_node
            for orig_id, group in self.group_original_to_node.items():
                if group.elementId == missing_group_id:
                    group_template = group
                    break

            if group_template:
                # Create a visual placeholder for this complex
                placeholder_group = copy.deepcopy(group_template)

                # Calculate position and size based on monomers
                monomer_xs = []
                monomer_ys = []
                for datanode in pathway_datanodes:
                    if hasattr(datanode, 'groupRef') and datanode.groupRef == missing_group_id:
                        if hasattr(datanode, 'graphics') and datanode.graphics:
                            monomer_xs.append(datanode.graphics.centerX)
                            monomer_ys.append(datanode.graphics.centerY)

                if monomer_xs and monomer_ys:
                    # Position group to overlay the monomers
                    center_x = sum(monomer_xs) / len(monomer_xs)
                    center_y = sum(monomer_ys) / len(monomer_ys)

                    # Calculate size
                    min_x = min(monomer_xs)
                    max_x = max(monomer_xs)
                    min_y = min(monomer_ys)
                    max_y = max(monomer_ys)

                    group_padding = 20
                    group_width = (max_x - min_x) + 100 + group_padding * 2
                    group_height = (max_y - min_y) + 60 + group_padding * 2

                    placeholder_group.graphics = Graphics(
                        centerX=center_x,
                        centerY=center_y,
                        width=group_width,
                        height=group_height,
                        textColor='666666',
                        fontName='Arial',
                        fontWeight=False,
                        fontStyle=False,
                        fontDecoration=False,
                        fontStrikethru=False,
                        fontSize=10.0,
                        hAlign=HAlign.CENTER,
                        vAlign=VAlign.MIDDLE,
                        borderColor='9900CC',
                        borderStyle=BorderStyle.DASHED,
                        borderWidth=2.0,
                        fillColor='F5E6FF',
                        shapeType=ShapeType.RECTANGLE,
                        zOrder=-1
                    )

                    pathway_groups.append(placeholder_group)

        return pathway_datanodes, pathway_groups, compound_node_map, placeholder_stats

    def _create_pathway_interactions(self, pathway_components, compound_node_map, pathway_datanodes):
        """
        Create all interactions for pathway.

        Args:
            pathway_components: Dictionary of pathway components
            compound_node_map: Map of compound IDs to positioned compound nodes
            pathway_datanodes: List of positioned DataNode objects in the pathway
        """
        pathway_interactions = []
        interaction_counter = 0

        # Create lookup map for positioned nodes (by elementId)
        positioned_node_map = {node.elementId: node for node in pathway_datanodes}

        # Extract primary info for easy access
        primary_info = pathway_components.get('primary_info', {})

        # Gene-protein interactions
        for gene_id in pathway_components['genes']:
            if gene_id in self.gene_protein_mapping:
                positioned_connections = []
                for connection in self.gene_protein_mapping[gene_id]:
                    gene_node = connection['gene_node']
                    protein_node = connection['protein_node']

                    positioned_gene = positioned_node_map.get(gene_node.elementId)
                    positioned_protein = positioned_node_map.get(protein_node.elementId)

                    if positioned_gene is not None and positioned_protein is not None:
                        positioned_connections.append({
                            'protein_id': connection['protein_id'],
                            'gene_node': positioned_gene,
                            'protein_node': positioned_protein
                        })

                if positioned_connections:
                    gene_protein_interactions = self._create_gene_protein_interactions(
                        gene_id, positioned_connections, interaction_counter
                    )
                    pathway_interactions.extend(gene_protein_interactions)
                    interaction_counter += 1

        # Protein-reaction interactions
        for reaction_data, original_reaction_id in pathway_components['reactions']:
            enzyme_nodes = []
            if original_reaction_id in self.reaction_to_enzymes:
                for enzyme_info in self.reaction_to_enzymes[original_reaction_id]:
                    protein_node = enzyme_info['protein_node']
                    if protein_node:
                        positioned_enzyme = positioned_node_map.get(protein_node.elementId, protein_node)
                        positioned_enzyme_info = enzyme_info.copy()
                        positioned_enzyme_info['protein_node'] = positioned_enzyme
                        enzyme_nodes.append(positioned_enzyme_info)
                    else:
                        enzyme_nodes.append(enzyme_info)

            # Get primary compound info for this reaction
            reaction_primary_info = primary_info.get(original_reaction_id, None)

            reaction_interactions = self._create_standard_reaction_with_central_anchor(
                reaction_data, compound_node_map, enzyme_nodes, original_reaction_id, positioned_node_map,
                primary_info=reaction_primary_info
            )
            pathway_interactions.extend(reaction_interactions)

        # Add regulation interactions (may create placeholder nodes and groups)
        regulation_interactions, new_placeholder_nodes, new_regulator_groups = self._create_regulation_interactions(
            pathway_components, compound_node_map, positioned_node_map
        )
        pathway_interactions.extend(regulation_interactions)

        # Add any new placeholder nodes and groups created for regulators
        pathway_datanodes.extend(new_placeholder_nodes)

        return pathway_interactions, new_regulator_groups

    def _create_regulation_interactions(self, pathway_components, compound_node_map, positioned_node_map):
        """
        Create regulation interactions for pathways.

        Handles cases where:
        1. Regulator is a compound -> draws from compound to reaction anchor
        2. Regulator is a protein -> draws from protein to reaction anchor
        3. Regulator is missing -> creates placeholder compound and draws regulation

        Args:
            pathway_components: Dictionary of pathway components
            compound_node_map: Map of compound IDs to positioned compound nodes
            positioned_node_map: Map of node element IDs to positioned nodes

        Returns:
            tuple: (list of regulation interactions, list of new placeholder nodes, list of new groups)
        """
        regulation_interactions = []
        new_placeholder_nodes = []
        new_regulator_groups = []  # Track new complex groups for regulators

        # Iterate through reactions in this pathway
        for reaction_data, original_reaction_id in pathway_components['reactions']:
            # Check if this reaction has regulations
            if original_reaction_id not in self.regulation_by_reaction:
                continue

            # Get the reaction's central anchor
            reaction_interaction = reaction_data['interaction']
            if not reaction_interaction.anchors:
                continue

            central_anchor = reaction_interaction.anchors[0]

            # Process each regulation for this reaction
            for reg_idx, reg_data in enumerate(self.regulation_by_reaction[original_reaction_id]):
                regulator_id = reg_data.get('regulator')
                if not regulator_id:
                    continue

                # Find the regulator node
                regulator_node = None

                # Try compound first (check in positioned nodes)
                if regulator_id in compound_node_map:
                    regulator_node = compound_node_map[regulator_id]
                # Try compound in loaded data
                elif regulator_id in self.compound_original_to_node:
                    # Compound exists in compounds.dat - add it as full DataNode
                    import copy
                    compound_entity = self.compound_original_to_node[regulator_id]
                    regulator_node = copy.deepcopy(compound_entity)

                    # Use standard metabolite graphics if missing
                    if not hasattr(regulator_node, 'graphics') or regulator_node.graphics is None:
                        regulator_node.graphics = standard_graphics.create_metabolite_graphics(0.0, 0.0)

                    # Position near the reaction (to the left)
                    if len(reaction_interaction.waypoints) >= 2:
                        reaction_center_x = (reaction_interaction.waypoints[0].x + reaction_interaction.waypoints[1].x) / 2
                        reaction_center_y = (reaction_interaction.waypoints[0].y + reaction_interaction.waypoints[1].y) / 2
                    else:
                        reaction_center_x = reaction_interaction.waypoints[0].x
                        reaction_center_y = reaction_interaction.waypoints[0].y

                    regulator_node.graphics.centerX = reaction_center_x - 150
                    regulator_node.graphics.centerY = reaction_center_y - 50

                    # Add to maps and track
                    compound_node_map[regulator_id] = regulator_node
                    positioned_node_map[regulator_node.elementId] = regulator_node
                    new_placeholder_nodes.append(regulator_node)
                # Try protein
                elif regulator_id in self.protein_original_to_node:
                    protein_entity = self.protein_original_to_node[regulator_id]
                    # Check if this protein isin this pathway
                    if protein_entity.elementId in positioned_node_map:
                        regulator_node = positioned_node_map[protein_entity.elementId]
                    else:
                        # Protein exists but not in this pathway - need to add it
                        import copy
                        from scripts.data_structure.wiki_data_structure import HAlign, VAlign, BorderStyle, ShapeType

                        # Check if this is a complex (Group object)
                        if isinstance(protein_entity, Group):
                            # It's a complex - we need to add the group AND its monomers
                            # For regulation line connection, we'll use one of the monomers
                            # Find monomers for this complex
                            if regulator_id in self.monomer_by_complex or protein_entity.elementId in self.monomer_by_complex:
                                monomer_key = regulator_id if regulator_id in self.monomer_by_complex else protein_entity.elementId
                                monomers = self.monomer_by_complex[monomer_key]

                                if monomers:
                                    if len(reaction_interaction.waypoints) >= 2:
                                        reaction_center_x = (reaction_interaction.waypoints[0].x + reaction_interaction.waypoints[1].x) / 2
                                        reaction_center_y = (reaction_interaction.waypoints[0].y + reaction_interaction.waypoints[1].y) / 2
                                    else:
                                        reaction_center_x = reaction_interaction.waypoints[0].x
                                        reaction_center_y = reaction_interaction.waypoints[0].y

                                    complex_x = reaction_center_x - 150
                                    complex_y = reaction_center_y - 50

                                    # Add monomers
                                    num_monomers = len(monomers)
                                    spacing = 90
                                    for i, monomer_template in enumerate(monomers):
                                        monomer_node = copy.deepcopy(monomer_template)

                                        if not hasattr(monomer_node, 'graphics') or monomer_node.graphics is None:
                                            monomer_node.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)

                                        offset_x = (i - (num_monomers - 1) / 2) * spacing
                                        monomer_node.graphics.centerX = complex_x + offset_x
                                        monomer_node.graphics.centerY = complex_y

                                        positioned_node_map[monomer_node.elementId] = monomer_node
                                        new_placeholder_nodes.append(monomer_node)

                                    # Use first monomer as regulator node for the regulation line
                                    first_monomer = copy.deepcopy(monomers[0])
                                    if not hasattr(first_monomer, 'graphics') or first_monomer.graphics is None:
                                        first_monomer.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)
                                    first_monomer.graphics.centerX = complex_x
                                    first_monomer.graphics.centerY = complex_y
                                    regulator_node = positioned_node_map.get(first_monomer.elementId, first_monomer)

                                    # Create and add the complex group overlay
                                    complex_group = copy.deepcopy(protein_entity)

                                    # Set group graphics to overlay all monomers
                                    total_width = spacing * (num_monomers - 1) + 100
                                    group_padding = 20
                                    group_width = total_width + group_padding * 2
                                    group_height = 60 + group_padding * 2

                                    complex_group.graphics = Graphics(
                                        centerX=complex_x,
                                        centerY=complex_y,
                                        width=group_width,
                                        height=group_height,
                                        textColor='666666',  # Gray text
                                        fontName='Arial',
                                        fontWeight=False,
                                        fontStyle=False,
                                        fontDecoration=False,
                                        fontStrikethru=False,
                                        fontSize=10.0,
                                        hAlign=HAlign.CENTER,
                                        vAlign=VAlign.MIDDLE,
                                        borderColor='9900CC',  # Purple border
                                        borderStyle=BorderStyle.DASHED,
                                        borderWidth=2.0,
                                        fillColor='F5E6FF',  # Very light purple fill
                                        shapeType=ShapeType.RECTANGLE,
                                        zOrder=-1  # Behind other elements
                                    )

                                    new_regulator_groups.append(complex_group)
                                else:
                                    # Complex has no monomers - skip
                                    continue
                            else:
                                # Complex not found in monomer mapping - skip
                                continue
                        else:
                            # Regular protein (not a complex)
                            regulator_node = copy.deepcopy(protein_entity)

                            # Use standard protein graphics if missing
                            if not hasattr(regulator_node, 'graphics') or regulator_node.graphics is None:
                                regulator_node.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)

                            # Position near the reaction
                            if len(reaction_interaction.waypoints) >= 2:
                                reaction_center_x = (reaction_interaction.waypoints[0].x + reaction_interaction.waypoints[1].x) / 2
                                reaction_center_y = (reaction_interaction.waypoints[0].y + reaction_interaction.waypoints[1].y) / 2
                            else:
                                reaction_center_x = reaction_interaction.waypoints[0].x
                                reaction_center_y = reaction_interaction.waypoints[0].y

                            regulator_node.graphics.centerX = reaction_center_x - 150
                            regulator_node.graphics.centerY = reaction_center_y - 50

                            # Add to maps and track
                            positioned_node_map[regulator_node.elementId] = regulator_node
                            new_placeholder_nodes.append(regulator_node)
                else:
                    # Regulator doesn't exist in compounds.dat or proteins.dat - skip this regulation
                    continue

                # Skip only if regulator node has no graphics (shouldn't happen now)
                if not regulator_node or not hasattr(regulator_node, 'graphics'):
                    continue

                # Get arrow head type from regulation mode
                mode = reg_data.get('mode')
                if mode == '+':
                    arrow_head = ArrowHeadType.STIMULATION
                elif mode == '-':
                    arrow_head = ArrowHeadType.INHIBITION
                else:
                    arrow_head = ArrowHeadType.DIRECTED

                # Calculate positions for regulation line
                # From regulator node (right side) to reaction center
                regulator_x = regulator_node.graphics.centerX + regulator_node.graphics.width / 2
                regulator_y = regulator_node.graphics.centerY

                # Calculate reaction center (midpoint between first two waypoints of reaction)
                if len(reaction_interaction.waypoints) >= 2:
                    reaction_center_x = (reaction_interaction.waypoints[0].x + reaction_interaction.waypoints[1].x) / 2
                    reaction_center_y = (reaction_interaction.waypoints[0].y + reaction_interaction.waypoints[1].y) / 2
                else:
                    reaction_center_x = reaction_interaction.waypoints[0].x
                    reaction_center_y = reaction_interaction.waypoints[0].y

                # Create regulation interaction
                regulation_interaction = Interaction(
                    elementId=self.id_manager.register_id(
                        f"{reaction_interaction.elementId}_regulation_{reg_idx}"
                    ),
                    waypoints=[
                        Point(
                            elementId=self.id_manager.register_id(
                                f"{reaction_interaction.elementId}_regulation_start_{reg_idx}"
                            ),
                            x=regulator_x,
                            y=regulator_y,
                            arrowHead=ArrowHeadType.UNDIRECTED,
                            elementRef=regulator_node.elementId,  # Connect to regulator DataNode
                            relX=1.0, relY=0.0
                        ),
                        Point(
                            elementId=self.id_manager.register_id(
                                f"{reaction_interaction.elementId}_regulation_end_{reg_idx}"
                            ),
                            x=reaction_center_x,
                            y=reaction_center_y,
                            arrowHead=arrow_head,
                            elementRef=central_anchor.elementId,  # Connect to reaction anchor
                            relX=0.0, relY=0.0
                        )
                    ],
                    anchors=[],
                    graphics=Graphics(
                        lineColor='000000',
                        lineStyle=LineStyle.SOLID,
                        lineWidth=1.0,
                        connectorType=ConnectorType.STRAIGHT
                    ),
                    comments=reg_data['interaction'].comments if 'interaction' in reg_data else [],
                    properties=reg_data['interaction'].properties if 'interaction' in reg_data else [],
                    citationRefs=reg_data['interaction'].citationRefs if 'interaction' in reg_data else []
                )

                regulation_interactions.append(regulation_interaction)

        return regulation_interactions, new_placeholder_nodes, new_regulator_groups

    def _expand_sub_pathways(self, reaction_list, visited=None):
        """
        Recursively expand sub-pathways in reaction list to get actual reactions.
        Returns a deduplicated list to avoid duplicate reaction IDs.

        Args:
            reaction_list (list): List of reaction IDs (which may include sub-pathway IDs)
            visited (set): Set of visited pathway IDs to prevent infinite recursion

        Returns:
            list: Expanded list with only unique actual reaction IDs (preserving order)
        """
        if visited is None:
            visited = set()

        expanded_reactions = []
        seen_reactions = set()  # Track seen reactions to avoid duplicates

        for item_id in reaction_list:
            # Check if this is a sub-pathway
            if item_id in self.pathway_records:
                # Prevent infinite recursion
                if item_id in visited:
                    continue
                visited.add(item_id)

                sub_record = self.pathway_records[item_id]
                sub_reaction_list = sub_record.get('REACTION-LIST', [])

                if isinstance(sub_reaction_list, str):
                    sub_reaction_list = [sub_reaction_list]
                elif not isinstance(sub_reaction_list, list):
                    sub_reaction_list = []

                if sub_reaction_list:
                    expanded = self._expand_sub_pathways(sub_reaction_list, visited)
                    for rxn_id in expanded:
                        if rxn_id not in seen_reactions:
                            expanded_reactions.append(rxn_id)
                            seen_reactions.add(rxn_id)
            else:
                if item_id not in seen_reactions:
                    expanded_reactions.append(item_id)
                    seen_reactions.add(item_id)

        return expanded_reactions

    def build_complete_pathway_with_genes(self, pathway_id, layout_type='grid', layout_params=None):
        """
        Build pathway with genes, proteins (including complexes), reactions, and citations.

        Args:
            pathway_id (str): ID of the pathway to build
            layout_type (str): Type of layout to use. Options: 'grid' (default) or 'forceatlas2'
            layout_params (dict, optional): Custom parameters for the layout algorithm
        """
        # Find pathway object
        pathway = None
        original_pathway_id = pathway_id

        for p in self.pathways:
            if (p.elementId == pathway_id or self.id_manager.get_sanitized_id(pathway_id) == p.elementId):
                pathway = p
                break

        if not pathway:
            raise ValueError(f"Pathway {pathway_id} not found")

        if original_pathway_id not in self.pathway_records:
            return pathway

        record = self.pathway_records[original_pathway_id]
        reaction_list = record.get('REACTION-LIST', [])

        if isinstance(reaction_list, str):
            reaction_list = [reaction_list]
        elif not isinstance(reaction_list, list):
            reaction_list = []

        if not reaction_list:
            return pathway

        expanded_reaction_list = self._expand_sub_pathways(reaction_list)

        if not expanded_reaction_list:
            return pathway

        pathway_components = self._collect_pathway_components(expanded_reaction_list, pathway_record=record)

        # Use selected layout algorithm
        if layout_type.lower() == 'forceatlas2':
            positions = self._calculate_component_positions_forceatlas2(pathway_components, [], layout_params)
        else:
            positions = self._calculate_component_positions(pathway_components)

        pathway_datanodes, pathway_groups, compound_node_map, placeholder_stats = self._create_pathway_datanodes(
            pathway_components, positions
        )


        pathway_interactions, new_regulator_groups = self._create_pathway_interactions(pathway_components, compound_node_map, pathway_datanodes)

        pathway.dataNodes = pathway_datanodes
        pathway.groups = pathway_groups + new_regulator_groups
        pathway.interactions = pathway_interactions

        # Collect element IDs for citation gathering
        element_ids = [original_pathway_id]

        for compound_id in pathway_components['compounds']:
            element_ids.append(compound_id)

        for protein_id in pathway_components['proteins']:
            element_ids.append(protein_id)

        for gene_id in pathway_components['genes']:
            element_ids.append(gene_id)

        for reaction_data, original_reaction_id in pathway_components['reactions']:
            element_ids.append(original_reaction_id)

            # Also add regulation IDs if this reaction has regulations
            if original_reaction_id in self.regulation_by_reaction:
                for reg_data in self.regulation_by_reaction[original_reaction_id]:
                    if 'interaction' in reg_data:
                        # Get the original regulation ID (before sanitization)
                        reg_interaction = reg_data['interaction']
                        for orig_id, sanit_id in self.id_manager.id_mapping.items():
                            if sanit_id == reg_interaction.elementId:
                                element_ids.append(orig_id)
                                break

        # Also collect citations from ALL datanodes actually in the pathway
        # This catches citations from placeholder nodes (regulators, etc.) that may have been added
        for datanode in pathway_datanodes:
            # Extract original ID from sanitized ID
            datanode_original_id = None
            for orig_id, sanit_id in self.id_manager.id_mapping.items():
                if sanit_id == datanode.elementId:
                    datanode_original_id = orig_id
                    break

            # Add if not already in element_ids
            if datanode_original_id and datanode_original_id not in element_ids:
                element_ids.append(datanode_original_id)

        # Also collect citations from ALL groups (complexes) in the pathway
        for group in pathway_groups:
            # Try to extract original ID from sanitized ID
            group_original_id = None
            for orig_id, sanit_id in self.id_manager.id_mapping.items():
                if sanit_id == group.elementId:
                    group_original_id = orig_id
                    break

            # If no mapping found, use the group's elementId directly
            if not group_original_id:
                group_original_id = group.elementId

            # Add if not already in element_ids
            if group_original_id and group_original_id not in element_ids:
                element_ids.append(group_original_id)

        pathway.citations = self.citation_manager.get_all_citations_for_pathway(element_ids)

        # Collect annotations for the pathway
        pathway_annotations = []
        referenced_annotation_ids = set()
        
        # Collect annotationRefs from all DataNodes
        for datanode in pathway_datanodes:
            if hasattr(datanode, 'annotationRefs') and datanode.annotationRefs:
                for ref in datanode.annotationRefs:
                    referenced_annotation_ids.add(ref.elementRef)
                    
        # Collect annotationRefs from Groups
        for group in pathway_groups:
            if hasattr(group, 'annotationRefs') and group.annotationRefs:
                for ref in group.annotationRefs:
                    referenced_annotation_ids.add(ref.elementRef)
                    
        # Add the actual Annotation objects
        for ann_id in referenced_annotation_ids:
            if ann_id in self.annotation_index:
                pathway_annotations.append(self.annotation_index[ann_id])
                
        pathway.annotations = pathway_annotations

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
            print(f"Found {len(missing_citation_refs)} missing citations, adding them...")
            from scripts.data_structure.wiki_data_structure import Citation, Xref
            for missing_ref in missing_citation_refs:
                # missing_ref is a sanitized elementId like "citation_PUB_12695547"
                # We need to reverse-engineer the original citation ID

                # Check if this citation already exists in citation_objects
                citation = None
                for orig_id, cit_obj in self.citation_manager.citation_objects.items():
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
                    citation = self.citation_manager.create_citation_object(unsanitized)

                # Add citation if we got one
                if citation and citation.elementId not in existing_citation_ids:
                    pathway.citations.append(citation)
                    existing_citation_ids.add(citation.elementId)

        self._update_pathway_board_size(pathway, pathway_datanodes)

        return pathway

    def deduplicate_pathway_elements(self, pathway):
        """
        Remove duplicate elements from a pathway by elementId.

        This ensures each element in the pathway has a unique elementId across ALL element types.
        When conflicts occur between different element types (e.g., Group and DataNode with same ID),
        priority is: Citations > DataNodes > Interactions > Groups

        Args:
            pathway: Pathway object to deduplicate

        Returns:
            Pathway: Deduplicated pathway
        """
        # Track all seen IDs across element types to detect cross-type conflicts
        all_seen_ids = set()
        # Track ID mappings for updating references
        id_remapping = {}

        # Step 1: Deduplicate Citations (highest priority)
        seen_citation_ids = set()
        deduplicated_citations = []
        for citation in pathway.citations:
            if citation.elementId not in seen_citation_ids:
                deduplicated_citations.append(citation)
                seen_citation_ids.add(citation.elementId)
                all_seen_ids.add(citation.elementId)
            else:
                print(f"WARNING: Removing duplicate Citation with elementId: {citation.elementId}")
        pathway.citations = deduplicated_citations

        # Step 2: Deduplicate DataNodes (check against all previous IDs)
        seen_datanode_ids = set()
        deduplicated_datanodes = []
        for datanode in pathway.dataNodes:
            if datanode.elementId not in seen_datanode_ids:
                if datanode.elementId in all_seen_ids:
                    # Cross-type conflict detected
                    print(f"WARNING: DataNode elementId '{datanode.elementId}' conflicts with existing element. Renaming...")
                    # Rename the conflicting datanode
                    original_id = datanode.elementId
                    datanode.elementId = self.id_manager.register_id(f"{original_id}_datanode")
                    id_remapping[original_id] = datanode.elementId
                    print(f"  Renamed DataNode '{original_id}' to '{datanode.elementId}'")
                deduplicated_datanodes.append(datanode)
                seen_datanode_ids.add(datanode.elementId)
                all_seen_ids.add(datanode.elementId)
            else:
                print(f"WARNING: Removing duplicate DataNode with elementId: {datanode.elementId}")
        pathway.dataNodes = deduplicated_datanodes

        # Step 3: Deduplicate Interactions (check against all previous IDs)
        seen_interaction_ids = set()
        deduplicated_interactions = []
        for interaction in pathway.interactions:
            if interaction.elementId not in seen_interaction_ids:
                if interaction.elementId in all_seen_ids:
                    # Cross-type conflict detected
                    print(f"WARNING: Interaction elementId '{interaction.elementId}' conflicts with existing element. Renaming...")
                    original_id = interaction.elementId
                    interaction.elementId = self.id_manager.register_id(f"{original_id}_interaction")
                    id_remapping[original_id] = interaction.elementId
                    print(f"  Renamed Interaction '{original_id}' to '{interaction.elementId}'")
                deduplicated_interactions.append(interaction)
                seen_interaction_ids.add(interaction.elementId)
                all_seen_ids.add(interaction.elementId)
            else:
                print(f"WARNING: Removing duplicate Interaction with elementId: {interaction.elementId}")
        pathway.interactions = deduplicated_interactions

        # Step 4: Deduplicate Groups (lowest priority, check against all previous IDs)
        seen_group_ids = set()
        deduplicated_groups = []
        for group in pathway.groups:
            if group.elementId not in seen_group_ids:
                if group.elementId in all_seen_ids:
                    # Cross-type conflict detected
                    print(f"WARNING: Group elementId '{group.elementId}' conflicts with existing element. Renaming...")
                    original_id = group.elementId
                    group.elementId = self.id_manager.register_id(f"{original_id}_group")
                    id_remapping[original_id] = group.elementId
                    print(f"  Renamed Group '{original_id}' to '{group.elementId}'")
                deduplicated_groups.append(group)
                seen_group_ids.add(group.elementId)
                all_seen_ids.add(group.elementId)
            else:
                print(f"WARNING: Removing duplicate Group with elementId: {group.elementId}")
        pathway.groups = deduplicated_groups

        # Step 5: Update references to renamed elements
        if id_remapping:
            print(f"\nUpdating {len(id_remapping)} element references...")

            # Update groupRef in DataNodes
            for datanode in pathway.dataNodes:
                if hasattr(datanode, 'groupRef') and datanode.groupRef and datanode.groupRef in id_remapping:
                    old_ref = datanode.groupRef
                    datanode.groupRef = id_remapping[old_ref]
                    print(f"  Updated DataNode '{datanode.elementId}' groupRef: {old_ref} -> {datanode.groupRef}")

            # Update groupRef in Interactions
            for interaction in pathway.interactions:
                if hasattr(interaction, 'groupRef') and interaction.groupRef and interaction.groupRef in id_remapping:
                    old_ref = interaction.groupRef
                    interaction.groupRef = id_remapping[old_ref]
                    print(f"  Updated Interaction '{interaction.elementId}' groupRef: {old_ref} -> {interaction.groupRef}")

            # Update groupRef in Groups (for nested groups)
            for group in pathway.groups:
                if hasattr(group, 'groupRef') and group.groupRef and group.groupRef in id_remapping:
                    old_ref = group.groupRef
                    group.groupRef = id_remapping[old_ref]
                    print(f"  Updated Group '{group.elementId}' groupRef: {old_ref} -> {group.groupRef}")

            # Update elementRef in Point waypoints
            for interaction in pathway.interactions:
                if hasattr(interaction, 'waypoints'):
                    for waypoint in interaction.waypoints:
                        if hasattr(waypoint, 'elementRef') and waypoint.elementRef and waypoint.elementRef in id_remapping:
                            old_ref = waypoint.elementRef
                            waypoint.elementRef = id_remapping[old_ref]
                            print(f"  Updated waypoint elementRef in Interaction '{interaction.elementId}': {old_ref} -> {waypoint.elementRef}")

        return pathway

    def _update_pathway_board_size(self, pathway, pathway_datanodes):
        """
        Update pathway board size to fit all components.
        """
        if pathway_datanodes:
            all_x = [node.graphics.centerX for node in pathway_datanodes if hasattr(node, 'graphics')]
            all_y = [node.graphics.centerY for node in pathway_datanodes if hasattr(node, 'graphics')]
            if all_x and all_y:
                pathway.graphics.boardWidth = max(1200, max(all_x) + 200)
                pathway.graphics.boardHeight = max(800, max(all_y) + 200)

    def find_all_pathways(self):
        """
        Find all pathways with reactions.
        Returns list of pathways with basic information.
        """
        pathways_list = []

        for pathway_id, pathway_record in self.pathway_records.items():
            reaction_list = pathway_record.get('REACTION-LIST', [])

            if isinstance(reaction_list, str):
                reaction_list = [reaction_list]
            elif not isinstance(reaction_list, list):
                reaction_list = []

            if not reaction_list:
                continue

            expanded_reaction_list = self._expand_sub_pathways(reaction_list)

            if not expanded_reaction_list:
                continue

            pathways_list.append({
                'pathway_id': pathway_id,
                'reactions': len(expanded_reaction_list),
                'title': pathway_record.get('COMMON-NAME', pathway_id)
            })

        return pathways_list

    def export_pathway_to_gpml(self, pathway, output_file):
        """
        Export pathway to GPML file with duplicate DataNode removal.
        """
        # Remove any duplicate DataNodes before exporting
        self._remove_duplicate_datanodes(pathway)

        writer = GPMLWriter()
        gpml_content = writer.write_pathway(pathway)
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(gpml_content)

    def _remove_duplicate_datanodes(self, pathway):
        """
        Remove duplicate DataNodes from a pathway (keeps first occurrence).
        This is a final safety check right before export to ensure no duplicates.

        Args:
            pathway: Pathway object to deduplicate
        """
        if not hasattr(pathway, 'dataNodes') or not pathway.dataNodes:
            return

        # Track seen IDs and remove duplicates, keeping first occurrence
        seen_ids = set()
        duplicates_removed = 0
        deduplicated = []

        for datanode in pathway.dataNodes:
            if datanode.elementId not in seen_ids:
                deduplicated.append(datanode)
                seen_ids.add(datanode.elementId)
            else:
                # Found a duplicate - skip it
                duplicates_removed += 1
                print(f"    [!] Removing duplicate DataNode: {datanode.elementId}")

        if duplicates_removed > 0:
            print(f"  [!] Removed {duplicates_removed} duplicate DataNode(s) from pathway")
            pathway.dataNodes = deduplicated

    def get_stats(self):
        """
        Get summary statistics about loaded data.
        """
        return {
            'compounds_loaded': len(self.compound_nodes),
            'genes_loaded': len(self.gene_nodes),
            'proteins_loaded': len(self.protein_nodes),
            'complexes_loaded': len(self.protein_groups),
            'reactions_loaded': len(self.reaction_data),
            'pathways_loaded': len(self.pathways),
            'gene_protein_mappings': len(self.gene_protein_mapping),
            'protein_reaction_mappings': len(self.protein_reaction_mapping),
            'reaction_enzyme_mappings': len(self.reaction_to_enzymes)
        }
