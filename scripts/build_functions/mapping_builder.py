
"""
Module for building biological entity mappings (Gene -> Protein -> Reaction).
"""
import os
from scripts.parsing_functions import parsing_utils

class MappingBuilder:
    def __init__(self, data_dir, genes_file, proteins_file, id_manager):
        self.data_dir = data_dir
        self.genes_file = genes_file
        self.proteins_file = proteins_file
        self.id_manager = id_manager
        
        # Mappings to be built
        self.gene_protein_mapping = {}
        self.protein_reaction_mapping = {}
        self.reaction_to_enzymes = {}
        self.regulation_by_reaction = {}
        self.regulation_by_regulator = {}
        self.regulation_by_regulated_entity = {}
        self.protein_to_genes = {}
        
        # Intermediate lookups
        self.enzrxn_to_reaction = {}
        self.monomer_by_complex = {}

    def build_all_mappings(self, gene_nodes_map, protein_nodes_map, protein_records, reaction_data, regulation_data):
        """
        Build all mappings and return them.
        """
        print("Building gene->protein mappings...")
        self._build_gene_protein_mapping(gene_nodes_map, protein_nodes_map)

        print("Building protein->reaction mappings...")
        self._build_complete_enzyme_mapping(protein_records, protein_nodes_map, reaction_data)

        print("Building regulation mappings...")
        self._build_regulation_mapping(regulation_data, protein_nodes_map) # compound map not needed for structure, just for stats maybe?

        self._build_protein_to_genes_mapping()
        print("Mappings complete!")
        
        return {
            'gene_protein_mapping': self.gene_protein_mapping,
            'protein_reaction_mapping': self.protein_reaction_mapping,
            'reaction_to_enzymes': self.reaction_to_enzymes,
            'regulation_by_reaction': self.regulation_by_reaction,
            'regulation_by_regulator': self.regulation_by_regulator,
            'regulation_by_regulated_entity': self.regulation_by_regulated_entity,
            'protein_to_genes': self.protein_to_genes,
            'monomer_by_complex': self.monomer_by_complex
        }

    def _build_gene_protein_mapping(self, gene_nodes_map, protein_nodes_map):
        """
        Build mapping from genes to their protein products.
        """
        # Load raw gene records
        genes_processor = parsing_utils.read_and_parse(self.genes_file)
        gene_by_id = {r.get('UNIQUE-ID'): r for r in genes_processor.records if r.get('UNIQUE-ID')}

        for gene_id, gene_record in gene_by_id.items():
            products = gene_record.get('PRODUCT', [])
            if not isinstance(products, list):
                products = [products] if products else []

            self.gene_protein_mapping[gene_id] = []

            for product_id in products:
                # Check if we have both the gene and protein nodes
                if gene_id in gene_nodes_map and product_id in protein_nodes_map:
                    self.gene_protein_mapping[gene_id].append({
                        'protein_id': product_id,
                        'gene_node': gene_nodes_map[gene_id],
                        'protein_node': protein_nodes_map[product_id]
                    })

    def _build_complete_enzyme_mapping(self, protein_records, protein_nodes_map, reaction_data):
        """
        Build complete protein->reaction mapping using enzrxns.dat bridge.
        """
        # Load enzrxns.dat
        enzrxns_file = os.path.join(self.data_dir, "enzrxns.dat")
        enzrxns_processor = parsing_utils.read_and_parse(enzrxns_file)

        for enzrxn_record in enzrxns_processor.records:
            enzrxn_id = enzrxn_record.get('UNIQUE-ID')
            reaction_id = enzrxn_record.get('REACTION')

            if enzrxn_id and reaction_id:
                if isinstance(reaction_id, list):
                    reaction_id = reaction_id[0] if reaction_id else None

                if reaction_id:
                    self.enzrxn_to_reaction[enzrxn_id] = reaction_id

        # Get available reaction IDs
        available_reaction_ids = set()
        for reaction_dict in reaction_data:
            reaction_id = reaction_dict['interaction'].elementId
            for orig_id, sanit_id in self.id_manager.id_mapping.items():
                if sanit_id == reaction_id:
                    available_reaction_ids.add(orig_id)
                    break

        # Build mapping
        for protein_id, protein_record in protein_records.items():
            if 'CATALYZES' in protein_record:
                types = protein_record.get('TYPES', [])
                if not isinstance(types, list): types = [types] if types else []
                is_complex = any('COMPLEX' in str(t).upper() for t in types)

                catalyzes = protein_record['CATALYZES']
                if not isinstance(catalyzes, list):
                    catalyzes = [catalyzes] if catalyzes else []

                protein_reactions = []

                for enzyme_id in catalyzes:
                    enzyme_id = str(enzyme_id).strip()

                    if enzyme_id in self.enzrxn_to_reaction:
                        target_reaction_id = self.enzrxn_to_reaction[enzyme_id]

                        if target_reaction_id in available_reaction_ids:
                            protein_reactions.append(target_reaction_id)

                            if target_reaction_id not in self.reaction_to_enzymes:
                                self.reaction_to_enzymes[target_reaction_id] = []

                            self.reaction_to_enzymes[target_reaction_id].append({
                                'protein_id': protein_id,
                                'protein_node': protein_nodes_map.get(protein_id),
                                'enzyme_id': enzyme_id,
                                'reaction_id': target_reaction_id,
                                'is_complex': is_complex
                            })

                if protein_reactions:
                    self.protein_reaction_mapping[protein_id] = protein_reactions

        # Store monomer nodes
        for node in protein_nodes_map.values():
            if hasattr(node, 'groupRef') and node.groupRef:
                complex_id = node.groupRef
                if complex_id not in self.monomer_by_complex:
                    self.monomer_by_complex[complex_id] = []
                self.monomer_by_complex[complex_id].append(node)

    def _build_regulation_mapping(self, regulation_data, protein_nodes_map):
        """
        Build mapping for regulation interactions.
        """
        # Load enzrxns if not already loaded (handled in constructor/previous method)
        
        for reg_data in regulation_data:
            regulated_entity = reg_data.get('regulated_entity')
            regulator = reg_data.get('regulator')

            if regulated_entity:
                if regulated_entity not in self.regulation_by_regulated_entity:
                    self.regulation_by_regulated_entity[regulated_entity] = []
                self.regulation_by_regulated_entity[regulated_entity].append(reg_data)

                if regulated_entity in self.enzrxn_to_reaction:
                    reaction_id = self.enzrxn_to_reaction[regulated_entity]
                    if reaction_id not in self.regulation_by_reaction:
                        self.regulation_by_reaction[reaction_id] = []
                    self.regulation_by_reaction[reaction_id].append(reg_data)

            if regulator:
                if regulator not in self.regulation_by_regulator:
                    self.regulation_by_regulator[regulator] = []
                self.regulation_by_regulator[regulator].append(reg_data)

    def _build_protein_to_genes_mapping(self):
        """
        Build a reverse mapping from proteins to genes.
        """
        for gene_id, connections in self.gene_protein_mapping.items():
            for connection in connections:
                protein_id = connection['protein_id']
                if protein_id not in self.protein_to_genes:
                    self.protein_to_genes[protein_id] = set()
                self.protein_to_genes[protein_id].add(gene_id)
