"""
Core Pathway Builder Module

This module provides the main CompletePathwayBuilderWithGenes class for loading
BioCyc data and building complete pathway networks with genes, proteins, complexes, and reactions.
"""

from scripts.data_structure.wiki_data_structure import (
    Pathway, DataNode, Interaction, Point, Graphics, ArrowHeadType, Anchor, AnchorShapeType,
    LineStyle, ConnectorType, Group, Comment, ShapeType, BorderStyle, DataNodeType,
    HAlign, VAlign
)
from scripts.build_functions.build_pathway_data_nodes import create_enhanced_pathways_from_file
from scripts.build_functions.build_compounds_data_nodes import create_enhanced_datanodes_from_compounds
from scripts.build_functions.build_protein_data_nodes import create_enhanced_datanodes_from_proteins
from scripts.build_functions.build_gene_data_nodes import create_enhanced_datanodes_from_genes
from scripts.build_functions.build_reaction_interaction import (
    create_all_reaction_interactions,
    create_standard_reaction_with_central_anchor
)
from scripts.build_functions.build_regulation_interaction import (
    create_all_regulation_interactions,
    create_pathway_regulation_interactions
)
from scripts.build_functions.citation_manager import CitationManager
from scripts.build_functions.mapping_builder import MappingBuilder
from scripts.object2gmpl.gpml_writer import GPMLWriter
from scripts.parsing_functions import parsing_utils
from scripts.utils import standard_graphics
from scripts.utils.layout import (
    calculate_component_positions,
    position_complex_elements,
    update_group_bounds
)
from scripts.utils.id_manager import IDManager
import math
import re
import copy
import os


class CompletePathwayBuilderWithGenes:
    """
    Pathway builder for BioCyc data with full gene-protein-reaction networks.

    This class loads all BioCyc data files and creates pathway visualizations
    showing the biological hierarchy: Genes -> Proteins/Complexes -> Reactions -> Compounds.

    """

    def __init__(self, compounds_file, genes_file, proteins_file, reactions_file, pathways_file, pubs_file="pubs.dat",
                 regulation_file="regulation.dat", version=None, organism_mapping=None):
        """
        Initialize builder by loading all BioCyc data files.
        """
        self.id_manager = IDManager()
        self.version = version
        self.organism_mapping = organism_mapping

        # Initialize CitationManager
        self.citation_manager = CitationManager(pubs_file)
        
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

        # Build biological mappings using MappingBuilder
        self._build_mappings()

    def _load_compounds(self, compounds_file):
        """Load and process compound data with citations."""
        self.compound_nodes = create_enhanced_datanodes_from_compounds(compounds_file, self.citation_manager)
        self._register_and_map_nodes(self.compound_nodes)
        self.compound_original_to_node = self._create_original_to_node_map(self.compound_nodes)

    def _load_genes(self, genes_file):
        """Load and process gene data with citations."""
        self.gene_nodes, self.gene_citations, self.gene_annotations = create_enhanced_datanodes_from_genes(
            genes_file, self.citation_manager, self.organism_mapping
        )
        # Index annotations
        for annotation in self.gene_annotations:
            self.annotation_index[annotation.elementId] = annotation
            
        self._register_and_map_nodes(self.gene_nodes)
        self.gene_original_to_node = self._create_original_to_node_map(self.gene_nodes)

    def _load_proteins(self, proteins_file):
        """Load and process protein data with complex support and enzyme information."""
        result = create_enhanced_datanodes_from_proteins(proteins_file, self.citation_manager, self.organism_mapping)

        if len(result) == 4:
            self.protein_nodes, self.protein_groups, self.protein_citations, self.protein_annotations = result
        elif len(result) == 3:
            self.protein_nodes, self.protein_groups, self.protein_citations = result
            self.protein_annotations = []
        else:
            self.protein_nodes, self.protein_citations = result
            self.protein_groups = []
            self.protein_annotations = []
            
        # Index annotations
        for annotation in self.protein_annotations:
            self.annotation_index[annotation.elementId] = annotation

        # Deduplicate
        seen_ids = set()
        unique_nodes = []
        for node in self.protein_nodes:
            if node.elementId not in seen_ids:
                unique_nodes.append(node)
                seen_ids.add(node.elementId)
        self.protein_nodes = unique_nodes

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

    def _load_regulations(self, regulation_file):
        """Load and process regulation data with citations."""
        self.regulation_data = create_all_regulation_interactions(regulation_file, self.citation_manager)
        self._register_regulation_ids()

    def _load_pathways(self, pathways_file):
        """Load and process pathway data with citations."""
        self.pathways = create_enhanced_pathways_from_file(pathways_file, self.citation_manager, self.version, self.organism_mapping)
        self._register_and_map_nodes(self.pathways)
        self.pathways_processor = parsing_utils.read_and_parse(pathways_file)
        self.pathway_records = {r['UNIQUE-ID']: r for r in self.pathways_processor.records if 'UNIQUE-ID' in r}

    def _build_mappings(self):
        """Build all biological relationship mappings using MappingBuilder."""
        mapping_builder = MappingBuilder(self.data_dir, self.genes_file, self.data_dir, self.id_manager)
        
        mappings = mapping_builder.build_all_mappings(
            self.gene_original_to_node,
            self.protein_original_to_node,
            self.protein_records,
            self.reaction_data,
            self.regulation_data
        )
        
        self.gene_protein_mapping = mappings['gene_protein_mapping']
        self.protein_reaction_mapping = mappings['protein_reaction_mapping']
        self.reaction_to_enzymes = mappings['reaction_to_enzymes']
        self.regulation_by_reaction = mappings['regulation_by_reaction']
        self.regulation_by_regulator = mappings['regulation_by_regulator']
        self.regulation_by_regulated_entity = mappings['regulation_by_regulated_entity']
        self.protein_to_genes = mappings['protein_to_genes']
        self.monomer_by_complex = mappings['monomer_by_complex']

    def _register_and_map_nodes(self, nodes):
        """Register and sanitize node IDs."""
        for node in nodes:
            original_id = node.elementId
            sanitized_id = self.id_manager.register_id(original_id)
            node.elementId = sanitized_id

    def _create_original_to_node_map(self, nodes):
        """Create mapping from original IDs to node objects."""
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

    def _create_gene_protein_interactions(self, gene_id, protein_connections, interaction_counter):
        """Create transcription/translation interactions from genes to proteins."""
        interactions = []
        for i, connection in enumerate(protein_connections):
            gene_node = connection['gene_node']
            protein_node = connection['protein_node']

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

            if not hasattr(gene_node, 'graphics') or gene_node.graphics is None:
                continue

            for p_idx, target_protein in enumerate(target_protein_nodes):
                if not hasattr(target_protein, 'graphics') or target_protein.graphics is None:
                    continue

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

    def _parse_primaries(self, primaries_field):
        """Parse PRIMARIES field from pathway record."""
        primary_map = {}
        if not primaries_field: return primary_map
        if not isinstance(primaries_field, list): primaries_field = [primaries_field]

        for entry in primaries_field:
            entry_str = str(entry)
            match = re.match(r'\("([^"]+)"\s+\(([^)]+)\)\s+\(([^)]+)\)\)', entry_str)
            if match:
                reaction_id = match.group(1)
                reactants_str = match.group(2)
                products_str = match.group(3)
                reactants = [c.strip().strip('"') for c in reactants_str.split() if c.strip().strip('"')]
                products = [c.strip().strip('"') for c in products_str.split() if c.strip().strip('"')]
                primary_map[reaction_id] = {'reactants': reactants, 'products': products}
        return primary_map

    def _parse_reaction_layout(self, layout_field):
        """Parse REACTION-LAYOUT field from pathway record."""
        layout_map = {}
        if not layout_field: return layout_map
        if not isinstance(layout_field, list): layout_field = [layout_field]

        for entry in layout_field:
            entry_str = str(entry)
            reaction_match = re.match(r'\(([^\s]+)', entry_str)
            if not reaction_match: continue
            reaction_id = reaction_match.group(1)

            left_match = re.search(r':LEFT-PRIMARIES\s+([^)]+)', entry_str)
            left_primaries = [c.strip() for c in left_match.group(1).split() if c.strip() and not c.startswith(':')] if left_match else []

            right_match = re.search(r':RIGHT-PRIMARIES\s+([^)]+)', entry_str)
            right_primaries = [c.strip() for c in right_match.group(1).split() if c.strip() and not c.startswith(':')] if right_match else []

            direction_match = re.search(r':DIRECTION\s+:([^\s)]+)', entry_str)
            direction = direction_match.group(1) if direction_match else 'L2R'

            layout_map[reaction_id] = {
                'left_primaries': left_primaries,
                'right_primaries': right_primaries,
                'direction': direction
            }
        return layout_map

    def _collect_pathway_components(self, reaction_list, pathway_record=None):
        """Collect all components (genes, proteins, compounds, reactions) for a pathway."""
        pathway_compounds = set()
        pathway_reactions = []
        pathway_proteins = set()
        pathway_genes = set()
        primary_info = {}

        if pathway_record:
            primaries_field = pathway_record.get('PRIMARIES', [])
            primaries_map = self._parse_primaries(primaries_field)
            layout_field = pathway_record.get('REACTION-LAYOUT', [])
            layout_map = self._parse_reaction_layout(layout_field)

            for rxn_id in set(list(primaries_map.keys()) + list(layout_map.keys())):
                primary_info[rxn_id] = {'left_primaries': [], 'right_primaries': []}
                if rxn_id in layout_map:
                    primary_info[rxn_id]['left_primaries'] = layout_map[rxn_id]['left_primaries']
                    primary_info[rxn_id]['right_primaries'] = layout_map[rxn_id]['right_primaries']
                    primary_info[rxn_id]['direction'] = layout_map[rxn_id].get('direction', 'L2R')
                elif rxn_id in primaries_map:
                    primary_info[rxn_id]['left_primaries'] = primaries_map[rxn_id]['reactants']
                    primary_info[rxn_id]['right_primaries'] = primaries_map[rxn_id]['products']

        for reaction_id in reaction_list:
            has_layout = reaction_id in primary_info and 'direction' in primary_info[reaction_id]
            if not has_layout:
                for path_id, path_record in self.pathway_records.items():
                    layout_field = path_record.get('REACTION-LAYOUT', [])
                    if layout_field:
                        layout_map = self._parse_reaction_layout(layout_field)
                        if reaction_id in layout_map:
                            primary_info[reaction_id] = {
                                'left_primaries': layout_map[reaction_id]['left_primaries'],
                                'right_primaries': layout_map[reaction_id]['right_primaries'],
                                'direction': layout_map[reaction_id].get('direction', 'L2R')
                            }
                            break
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
                                break

        for original_reaction_id in reaction_list:
            sanitized_reaction_id = self.id_manager.get_sanitized_id(original_reaction_id)
            if sanitized_reaction_id in self.reaction_index:
                reaction_data = self.reaction_index[sanitized_reaction_id]
                pathway_reactions.append((reaction_data, original_reaction_id))
                for reactant in reaction_data['reactants']: pathway_compounds.add(reactant['compound_id'])
                for product in reaction_data['products']: pathway_compounds.add(product['compound_id'])
                if original_reaction_id in self.reaction_to_enzymes:
                    for enzyme_info in self.reaction_to_enzymes[original_reaction_id]:
                        pathway_proteins.add(enzyme_info['protein_id'])

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

    def _create_pathway_datanodes(self, pathway_components, positions):
        """Create DataNode objects for pathway components."""
        pathway_datanodes = []
        pathway_groups = []
        compound_node_map = {}
        placeholder_stats = {'compounds': []}

        # Genes
        for gene_id in pathway_components['genes']:
            if gene_id in self.gene_original_to_node:
                gene_node = copy.deepcopy(self.gene_original_to_node[gene_id])
                if not hasattr(gene_node, 'graphics') or gene_node.graphics is None:
                    gene_node.graphics = standard_graphics.create_gene_graphics(0.0, 0.0)
                if gene_id in positions['genes']:
                    x, y = positions['genes'][gene_id]
                    gene_node.graphics.centerX = x
                    gene_node.graphics.centerY = y
                pathway_datanodes.append(gene_node)

        # Proteins
        processed_proteins = set()
        global_monomer_ids = set()
        for protein_id in pathway_components['proteins']:
            protein_entity = self.protein_original_to_node.get(protein_id)
            is_complex_group = protein_id in self.group_original_to_node or isinstance(protein_entity, Group)

            if is_complex_group:
                if protein_id in self.group_original_to_node:
                    complex_group = copy.deepcopy(self.group_original_to_node[protein_id])
                else:
                    complex_group = copy.deepcopy(protein_entity)

                complex_x, complex_y = positions['proteins'].get(protein_id, (150, 150))
                monomer_nodes_in_complex = []
                complex_monomer_ids = set()

                for node in self.protein_nodes:
                    if hasattr(node, 'groupRef') and (node.groupRef == complex_group.elementId or 
                       (protein_id in self.id_manager.id_mapping and node.groupRef == self.id_manager.id_mapping[protein_id])):
                        if node.elementId in global_monomer_ids or node.elementId in complex_monomer_ids: continue
                        
                        monomer_node = copy.deepcopy(node)
                        if not hasattr(monomer_node, 'graphics') or monomer_node.graphics is None:
                            monomer_node.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)
                        
                        monomer_nodes_in_complex.append(monomer_node)
                        complex_monomer_ids.add(node.elementId)
                        global_monomer_ids.add(node.elementId)
                        processed_proteins.add(node.elementId)

                position_complex_elements(complex_group, monomer_nodes_in_complex, complex_x, complex_y)
                
                # Set static graphics properties
                complex_group.graphics.textColor = '666666'
                complex_group.graphics.fontName = 'Arial'
                complex_group.graphics.hAlign = HAlign.CENTER
                complex_group.graphics.vAlign = VAlign.MIDDLE
                complex_group.graphics.borderColor = '9900CC'
                complex_group.graphics.borderStyle = BorderStyle.DASHED
                complex_group.graphics.borderWidth = 2.0
                complex_group.graphics.fillColor = 'F5E6FF'
                complex_group.graphics.shapeType = ShapeType.RECTANGLE
                complex_group.graphics.zOrder = -1

                pathway_datanodes.extend(monomer_nodes_in_complex)
                pathway_groups.append(complex_group)

            elif protein_id in self.protein_original_to_node and protein_id not in processed_proteins:
                protein_entity = self.protein_original_to_node[protein_id]
                if not isinstance(protein_entity, Group):
                    protein_node = copy.deepcopy(protein_entity)
                    if not hasattr(protein_node, 'graphics') or protein_node.graphics is None:
                        protein_node.graphics = standard_graphics.create_protein_graphics(0.0, 0.0)
                    
                    x, y = positions['proteins'].get(protein_id, (150, 150))
                    protein_node.graphics.centerX = x
                    protein_node.graphics.centerY = y
                    pathway_datanodes.append(protein_node)
                    processed_proteins.add(protein_id)

        # Compounds
        for compound_id in pathway_components['compounds']:
            if compound_id in self.compound_original_to_node:
                compound_node = copy.deepcopy(self.compound_original_to_node[compound_id])
                if not hasattr(compound_node, 'graphics') or compound_node.graphics is None:
                    compound_node.graphics = standard_graphics.create_metabolite_graphics(0.0, 0.0)
                
                x, y = positions['compounds'].get(compound_id, (200, 250))
                compound_node.graphics.centerX = x
                compound_node.graphics.centerY = y
                pathway_datanodes.append(compound_node)
                compound_node_map[compound_id] = compound_node
            else:
                from scripts.data_structure.wiki_data_structure import DataNode, Xref, DataNodeType, Property
                from scripts.utils.HTML_cleaner import clean_text_label
                
                sanitized_id = self.id_manager.register_id(compound_id)
                placeholder_node = DataNode(
                    elementId=sanitized_id, textLabel=clean_text_label(compound_id),
                    type=DataNodeType.METABOLITE, xref=Xref(identifier=compound_id, dataSource=''),
                    graphics=standard_graphics.create_metabolite_graphics(0.0, 0.0),
                    properties=[Property(key='IsPlaceholder', value='true')]
                )
                placeholder_node.graphics.fillColor = 'F0F0F0'
                placeholder_node.graphics.borderColor = 'AAAAAA'
                
                x, y = positions['compounds'].get(compound_id, (200, 250))
                placeholder_node.graphics.centerX = x
                placeholder_node.graphics.centerY = y
                pathway_datanodes.append(placeholder_node)
                compound_node_map[compound_id] = placeholder_node
                placeholder_stats['compounds'].append(compound_id)

        # Handle missing groups
        group_ids_in_pathway = {group.elementId for group in pathway_groups}
        datanode_ids_in_pathway = {node.elementId for node in pathway_datanodes}
        missing_group_refs = set()
        for datanode in pathway_datanodes:
            if hasattr(datanode, 'groupRef') and datanode.groupRef and datanode.groupRef not in group_ids_in_pathway:
                missing_group_refs.add(datanode.groupRef)

        nested_complexes = missing_group_refs & datanode_ids_in_pathway
        for datanode in pathway_datanodes:
            if hasattr(datanode, 'groupRef') and datanode.groupRef in nested_complexes:
                datanode.groupRef = None
        
        missing_group_refs = missing_group_refs - nested_complexes
        
        for missing_group_id in missing_group_refs:
            group_template = None
            for orig_id, group in self.group_original_to_node.items():
                if group.elementId == missing_group_id:
                    group_template = group
                    break
            
            if group_template:
                placeholder_group = copy.deepcopy(group_template)
                members = [n for n in pathway_datanodes if hasattr(n, 'groupRef') and n.groupRef == missing_group_id]
                update_group_bounds(placeholder_group, members)
                
                if members:
                    placeholder_group.graphics.textColor = '666666'
                    placeholder_group.graphics.fontName = 'Arial'
                    placeholder_group.graphics.hAlign = HAlign.CENTER
                    placeholder_group.graphics.vAlign = VAlign.MIDDLE
                    placeholder_group.graphics.borderColor = '9900CC'
                    placeholder_group.graphics.borderStyle = BorderStyle.DASHED
                    placeholder_group.graphics.borderWidth = 2.0
                    placeholder_group.graphics.fillColor = 'F5E6FF'
                    placeholder_group.graphics.shapeType = ShapeType.RECTANGLE
                    placeholder_group.graphics.zOrder = -1
                    pathway_groups.append(placeholder_group)

        return pathway_datanodes, pathway_groups, compound_node_map, placeholder_stats

    def _create_pathway_interactions(self, pathway_components, compound_node_map, pathway_datanodes):
        """Create all interactions for pathway."""
        pathway_interactions = []
        interaction_counter = 0
        positioned_node_map = {node.elementId: node for node in pathway_datanodes}
        primary_info = pathway_components.get('primary_info', {})

        # Gene-protein interactions
        for gene_id in pathway_components['genes']:
            if gene_id in self.gene_protein_mapping:
                positioned_connections = []
                for connection in self.gene_protein_mapping[gene_id]:
                    positioned_gene = positioned_node_map.get(connection['gene_node'].elementId)
                    positioned_protein = positioned_node_map.get(connection['protein_node'].elementId)
                    if positioned_gene and positioned_protein:
                        positioned_connections.append({
                            'protein_id': connection['protein_id'],
                            'gene_node': positioned_gene,
                            'protein_node': positioned_protein
                        })
                
                if positioned_connections:
                    interactions = self._create_gene_protein_interactions(gene_id, positioned_connections, interaction_counter)
                    pathway_interactions.extend(interactions)
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

            reaction_interactions = create_standard_reaction_with_central_anchor(
                reaction_data, compound_node_map, enzyme_nodes, original_reaction_id, 
                positioned_node_map, self.id_manager, primary_info.get(original_reaction_id)
            )
            pathway_interactions.extend(reaction_interactions)

        # Regulation interactions
        regulation_interactions, new_nodes, new_groups = create_pathway_regulation_interactions(
            pathway_components, self.regulation_by_reaction, compound_node_map, positioned_node_map,
            self.compound_original_to_node, self.protein_original_to_node, self.monomer_by_complex, self.id_manager
        )
        pathway_interactions.extend(regulation_interactions)
        pathway_datanodes.extend(new_nodes)

        return pathway_interactions, new_groups

    def _expand_sub_pathways(self, reaction_list, visited=None):
        """Recursively expand sub-pathways."""
        if visited is None: visited = set()
        expanded_reactions = []
        seen_reactions = set()

        for item_id in reaction_list:
            if item_id in self.pathway_records:
                if item_id in visited: continue
                visited.add(item_id)
                
                sub_reactions = self.pathway_records[item_id].get('REACTION-LIST', [])
                if isinstance(sub_reactions, str): sub_reactions = [sub_reactions]
                
                if sub_reactions:
                    for rxn_id in self._expand_sub_pathways(sub_reactions, visited):
                        if rxn_id not in seen_reactions:
                            expanded_reactions.append(rxn_id)
                            seen_reactions.add(rxn_id)
            else:
                if item_id not in seen_reactions:
                    expanded_reactions.append(item_id)
                    seen_reactions.add(item_id)
        return expanded_reactions

    def build_complete_pathway_with_genes(self, pathway_id):
        """Build pathway with genes, proteins, complexes, reactions, and citations."""
        pathway = None
        original_pathway_id = pathway_id
        for p in self.pathways:
            if (p.elementId == pathway_id or self.id_manager.get_sanitized_id(pathway_id) == p.elementId):
                pathway = p
                break
        
        if not pathway or original_pathway_id not in self.pathway_records:
            return pathway

        record = self.pathway_records[original_pathway_id]
        reaction_list = record.get('REACTION-LIST', [])
        if isinstance(reaction_list, str): reaction_list = [reaction_list]
        
        expanded_reaction_list = self._expand_sub_pathways(reaction_list) if reaction_list else []
        if not expanded_reaction_list: return pathway

        pathway_components = self._collect_pathway_components(expanded_reaction_list, pathway_record=record)
        positions = calculate_component_positions(pathway_components)
        
        pathway_datanodes, pathway_groups, compound_node_map, _ = self._create_pathway_datanodes(pathway_components, positions)
        pathway_interactions, new_regulator_groups = self._create_pathway_interactions(pathway_components, compound_node_map, pathway_datanodes)

        pathway.dataNodes = pathway_datanodes
        pathway.groups = pathway_groups + new_regulator_groups
        pathway.interactions = pathway_interactions

        # Citations
        element_ids = [original_pathway_id]
        element_ids.extend(pathway_components['compounds'])
        element_ids.extend(pathway_components['proteins'])
        element_ids.extend(pathway_components['genes'])
        
        for _, rxn_id in pathway_components['reactions']:
            element_ids.append(rxn_id)
            if rxn_id in self.regulation_by_reaction:
                for reg in self.regulation_by_reaction[rxn_id]:
                    if 'interaction' in reg:
                        for orig, sanit in self.id_manager.id_mapping.items():
                            if sanit == reg['interaction'].elementId:
                                element_ids.append(orig)
                                break
        
        # Add original IDs from sanitized DataNodes
        for node in pathway_datanodes:
            for orig, sanit in self.id_manager.id_mapping.items():
                if sanit == node.elementId and orig not in element_ids:
                    element_ids.append(orig)
                    break
                    
        # Add original IDs from groups
        for group in pathway_groups:
            orig_id = None
            for orig, sanit in self.id_manager.id_mapping.items():
                if sanit == group.elementId:
                    orig_id = orig
                    break
            if not orig_id: orig_id = group.elementId
            if orig_id not in element_ids: element_ids.append(orig_id)

        pathway.citations = self.citation_manager.get_all_citations_for_pathway(element_ids)

        # Annotations
        pathway_annotations = []
        referenced_anns = set()
        for obj in pathway_datanodes + pathway_groups:
            if hasattr(obj, 'annotationRefs') and obj.annotationRefs:
                for ref in obj.annotationRefs: referenced_anns.add(ref.elementRef)
        
        for ann_id in referenced_anns:
            if ann_id in self.annotation_index:
                pathway_annotations.append(self.annotation_index[ann_id])
        pathway.annotations = pathway_annotations

        # Check for missing citations
        cited_refs = set()
        for obj in pathway_datanodes + pathway_groups + pathway_interactions:
            if hasattr(obj, 'citationRefs') and obj.citationRefs:
                for ref in obj.citationRefs: cited_refs.add(ref.elementRef)
        
        existing_citations = {c.elementId for c in pathway.citations if c.elementId}
        missing_refs = cited_refs - existing_citations
        
        if missing_refs:
            print(f"Found {len(missing_refs)} missing citations, adding them...")
            from scripts.data_structure.wiki_data_structure import Citation, Xref
            for missing in missing_refs:
                citation = None
                for _, cit_obj in self.citation_manager.citation_objects.items():
                    if cit_obj.elementId == missing:
                        citation = cit_obj
                        break
                
                if not citation:
                    unsanitized = missing.replace('citation_', '', 1)
                    if unsanitized.startswith('PUB_'): unsanitized = unsanitized.replace('_', '-', 1)
                    elif unsanitized.startswith('cit_'): unsanitized = unsanitized.replace('cit_', '', 1)
                    citation = self.citation_manager.create_citation_object(unsanitized)
                
                if citation and citation.elementId not in existing_citations:
                    pathway.citations.append(citation)
                    existing_citations.add(citation.elementId)

        self._update_pathway_board_size(pathway, pathway_datanodes)
        return pathway

    def _update_pathway_board_size(self, pathway, pathway_datanodes):
        """Update pathway board size to fit all components."""
        if pathway_datanodes:
            all_x = [node.graphics.centerX for node in pathway_datanodes if hasattr(node, 'graphics')]
            all_y = [node.graphics.centerY for node in pathway_datanodes if hasattr(node, 'graphics')]
            if all_x and all_y:
                pathway.graphics.boardWidth = max(1200, max(all_x) + 200)
                pathway.graphics.boardHeight = max(800, max(all_y) + 200)

    def deduplicate_pathway_elements(self, pathway):
        """Remove duplicate elements from a pathway by elementId."""
        all_ids = set()
        id_remap = {}

        # Citations
        unique_citations = []
        for cit in pathway.citations:
            if cit.elementId not in all_ids:
                unique_citations.append(cit)
                all_ids.add(cit.elementId)
        pathway.citations = unique_citations

        # DataNodes
        unique_datanodes = []
        for node in pathway.dataNodes:
            if node.elementId in all_ids:
                orig = node.elementId
                node.elementId = self.id_manager.register_id(f"{orig}_datanode")
                id_remap[orig] = node.elementId
            unique_datanodes.append(node)
            all_ids.add(node.elementId)
        pathway.dataNodes = unique_datanodes

        # Interactions
        unique_interactions = []
        for intr in pathway.interactions:
            if intr.elementId in all_ids:
                orig = intr.elementId
                intr.elementId = self.id_manager.register_id(f"{orig}_interaction")
                id_remap[orig] = intr.elementId
            unique_interactions.append(intr)
            all_ids.add(intr.elementId)
        pathway.interactions = unique_interactions

        # Groups
        unique_groups = []
        for grp in pathway.groups:
            if grp.elementId in all_ids:
                orig = grp.elementId
                grp.elementId = self.id_manager.register_id(f"{orig}_group")
                id_remap[orig] = grp.elementId
            unique_groups.append(grp)
            all_ids.add(grp.elementId)
        pathway.groups = unique_groups

        # Remap references
        if id_remap:
            for obj in pathway.dataNodes + pathway.interactions + pathway.groups:
                if hasattr(obj, 'groupRef') and obj.groupRef in id_remap:
                    obj.groupRef = id_remap[obj.groupRef]
            
            for intr in pathway.interactions:
                if hasattr(intr, 'waypoints'):
                    for wp in intr.waypoints:
                        if hasattr(wp, 'elementRef') and wp.elementRef in id_remap:
                            wp.elementRef = id_remap[wp.elementRef]

        return pathway

    def export_pathway_to_gpml(self, pathway, output_file):
        """Export pathway to GPML file."""
        self._remove_duplicate_datanodes(pathway)
        writer = GPMLWriter()
        gpml_content = writer.write_pathway(pathway)
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(gpml_content)

    def _remove_duplicate_datanodes(self, pathway):
        """Remove duplicate DataNodes from a pathway (keeps first occurrence)."""
        if not hasattr(pathway, 'dataNodes') or not pathway.dataNodes:
            return

        seen_ids = set()
        deduplicated = []

        for datanode in pathway.dataNodes:
            if datanode.elementId not in seen_ids:
                deduplicated.append(datanode)
                seen_ids.add(datanode.elementId)
        pathway.dataNodes = deduplicated

    def get_stats(self):
        """Get summary statistics."""
        return {
            'compounds_loaded': len(self.compound_nodes),
            'genes_loaded': len(self.gene_nodes),
            'proteins_loaded': len(self.protein_nodes),
            'reactions_loaded': len(self.reaction_data),
            'pathways_loaded': len(self.pathways),
        }

    def find_all_pathways(self):
        """Find all pathways with reactions."""
        pathways_list = []
        for pathway_id, pathway_record in self.pathway_records.items():
            reaction_list = pathway_record.get('REACTION-LIST', [])
            if isinstance(reaction_list, str): reaction_list = [reaction_list]
            if not isinstance(reaction_list, list): reaction_list = []
            
            if not reaction_list: continue
            
            expanded_reaction_list = self._expand_sub_pathways(reaction_list)
            if not expanded_reaction_list: continue
            
            pathways_list.append({
                'pathway_id': pathway_id,
                'reactions': len(expanded_reaction_list),
                'title': pathway_record.get('COMMON-NAME', pathway_id)
            })
        return pathways_list
