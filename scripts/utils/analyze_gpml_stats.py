#!/usr/bin/env python3
"""
GPML Statistics Analyzer
=========================

Analyzes all generated GPML files (pathways and reactions) to extract statistics:
- Interaction types (gene→protein, protein→reaction, etc.)
- DataNode types and their Xrefs
- Placeholders (total and unique)
- Publications (with/without PubMed ID or DOI)
- Coverage statistics (unique vs total counts)

Usage:
    python analyze_gpml_stats.py <output_dir>

Example:
    python analyze_gpml_stats.py ./output/biocyc_pathways_20251022_103448
"""

import sys
import os
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter
from datetime import datetime


class GPMLAnalyzer:
    """Analyzes GPML files for statistics."""

    def __init__(self):
        # Interaction statistics
        self.interaction_types = defaultdict(int)  # Type -> count
        self.interaction_type_xrefs = defaultdict(lambda: defaultdict(int))  # Type -> {dataSource -> count}
        self.unique_interactions = set()  # Unique interaction IDs

        # DataNode statistics
        self.datanode_types = defaultdict(int)  # Type -> count
        self.datanode_type_xrefs = defaultdict(lambda: defaultdict(int))  # Type -> {dataSource -> count}
        self.datanode_type_xrefs_unique = defaultdict(lambda: defaultdict(set))  # Type -> {dataSource -> set of unique IDs}
        self.unique_datanodes = set()  # Unique datanode IDs (textLabel)

        # Complex statistics (Groups that are complexes)
        self.total_complexes = 0
        self.unique_complexes = set()  # Unique complex IDs

        # Placeholder statistics
        self.total_placeholders = 0
        self.unique_placeholders = set()  # Unique placeholder IDs

        # Publication statistics
        self.total_citations = 0
        self.unique_citations = set()  # Unique citation IDs
        self.citations_with_pubmed = 0
        self.unique_citations_with_pubmed = set()
        self.citations_with_doi = 0
        self.unique_citations_with_doi = set()
        self.citations_without_either = 0
        self.unique_citations_without_either = set()

        # Taxonomy annotation statistics
        self.pathway_organisms = set()            # Unique pathway-level organism values
        self.pathway_organism_counts = defaultdict(int)  # organism -> file count
        # Per DataNode type: species name -> occurrence count
        self.taxonomy_by_type = defaultdict(lambda: defaultdict(int))
        # Per DataNode type: total nodes vs nodes with a taxonomy annotation
        self.total_by_type = defaultdict(int)
        self.annotated_by_type = defaultdict(int)

        # Pathway vs Reaction file tracking
        self.pathway_files = 0
        self.reaction_files = 0

        # Separate stats for pathways vs reactions
        self.pathway_stats = self._create_stats_dict()
        self.reaction_stats = self._create_stats_dict()

        # Keep legacy field so existing callers don't break
        self.protein_species = set()
        self.protein_species_counts = defaultdict(int)

        # File tracking
        self.files_processed = 0
        self.files_failed = 0

        # Property statistics
        self.property_keys = defaultdict(int)  # Property key -> count
        self.datanode_property_keys = defaultdict(int)  # DataNode property keys
        self.interaction_property_keys = defaultdict(int)  # Interaction property keys
        self.total_properties = 0

    def _create_stats_dict(self):
        """Create a statistics dictionary for pathways or reactions."""
        return {
            'datanode_types': defaultdict(int),
            'datanode_type_xrefs': defaultdict(lambda: defaultdict(int)),
            'unique_datanodes': set(),
            'interaction_types': defaultdict(int),
            'protein_species': set(),
            'protein_species_counts': defaultdict(int),
        }

    def parse_gpml_file(self, filepath):
        """Parse a single GPML file and extract statistics."""
        try:
            tree = ET.parse(filepath)
            root = tree.getroot()

            # Extract namespace
            ns = {'gpml': 'http://pathvisio.org/GPML/2021'}

            # Determine if pathway or reaction file
            is_reaction = 'individual_reactions' in filepath or 'single_reactions' in filepath
            file_stats = self.reaction_stats if is_reaction else self.pathway_stats

            if is_reaction:
                self.reaction_files += 1
            else:
                self.pathway_files += 1

            # Build per-file annotation index: elementId -> {value, type}
            annotation_index = self._build_annotation_index(root, ns)

            # Process Pathway organism
            self._process_pathway_organism(root, ns)

            # Process DataNodes
            self._process_datanodes(root, ns, file_stats, annotation_index)

            # Process Groups (for complexes)
            self._process_groups(root, ns, file_stats)

            # Process Interactions
            self._process_interactions(root, ns, file_stats)

            # Process Citations
            self._process_citations(root, ns)

            # Process Properties
            self._process_properties(root, ns)

            self.files_processed += 1

        except Exception as e:
            self.files_failed += 1
            print(f"  Error processing {os.path.basename(filepath)}: {e}")

    def _build_annotation_index(self, root, ns):
        """Build a per-file index of Annotation elements: elementId -> {value, type}."""
        index = {}
        for ann in root.findall('.//gpml:Annotation', ns):
            ann_id = ann.get('elementId', '')
            if ann_id:
                index[ann_id] = {
                    'value': ann.get('value', ''),
                    'type': ann.get('type', ''),
                }
        return index

    def _process_pathway_organism(self, root, ns):
        """Process Pathway organism attribute."""
        organism = root.get('organism', None)
        if organism:
            self.pathway_organisms.add(organism.strip())
            self.pathway_organism_counts[organism.strip()] += 1

    def _process_datanodes(self, root, ns, file_stats, annotation_index=None):
        """Process DataNode elements."""
        if annotation_index is None:
            annotation_index = {}

        for datanode in root.findall('.//gpml:DataNode', ns):
            node_type = datanode.get('type', 'Unknown')
            self.datanode_types[node_type] += 1
            file_stats['datanode_types'][node_type] += 1
            self.total_by_type[node_type] += 1

            text_label = datanode.get('textLabel', '')
            element_id = datanode.get('elementId', '')
            unique_id = text_label if text_label else element_id
            if unique_id:
                self.unique_datanodes.add((node_type, unique_id))
                file_stats['unique_datanodes'].add((node_type, unique_id))

            # Check for placeholder
            for prop in datanode.findall('.//gpml:Property', ns):
                if prop.get('key') == 'IsPlaceholder' and prop.get('value') == 'true':
                    self.total_placeholders += 1
                    self.unique_placeholders.add((node_type, unique_id))
                    break

            # Get Xref
            xref = datanode.find('.//gpml:Xref', ns)
            if xref is not None:
                data_source = xref.get('dataSource', 'Unknown')
                if data_source and data_source != 'Unknown':
                    self.datanode_type_xrefs[node_type][data_source] += 1
                    file_stats['datanode_type_xrefs'][node_type][data_source] += 1
                    if unique_id:
                        self.datanode_type_xrefs_unique[node_type][data_source].add(unique_id)

            # Taxonomy annotations via AnnotationRef -> Annotation lookup
            has_taxonomy = False
            for ann_ref in datanode.findall('gpml:AnnotationRef', ns):
                ref_id = ann_ref.get('elementRef', '')
                ann = annotation_index.get(ref_id, {})
                if ann.get('type', '').lower() == 'taxonomy':
                    species_name = ann.get('value', '')
                    if species_name:
                        self.taxonomy_by_type[node_type][species_name] += 1
                        has_taxonomy = True
            if has_taxonomy:
                self.annotated_by_type[node_type] += 1

    def _process_groups(self, root, ns, file_stats):
        """Process Group elements (protein complexes)."""
        for group in root.findall('.//gpml:Group', ns):
            # Check if it's a complex
            is_complex = False
            group_type = group.get('type', '')

            # Check type attribute
            if group_type == 'Complex':
                is_complex = True

            # Check IsComplex property
            if not is_complex:
                for prop in group.findall('.//gpml:Property', ns):
                    if prop.get('key') == 'IsComplex' and prop.get('value') == 'true':
                        is_complex = True
                        break

            if is_complex:
                self.total_complexes += 1

                # Get unique identifier
                element_id = group.get('elementId', '')
                text_label = group.get('textLabel', '')
                unique_id = text_label if text_label else element_id

                if unique_id:
                    self.unique_complexes.add(unique_id)

    def _process_interactions(self, root, ns, file_stats):
        """Process Interaction elements."""
        for interaction in root.findall('.//gpml:Interaction', ns):
            element_id = interaction.get('elementId', '')

            # Skip anchor interactions (compound-to-anchor connections)
            if '_anchor' in element_id:
                # Check if this connects TO an anchor point
                skip_interaction = False
                for point in interaction.findall('.//gpml:Point', ns):
                    element_ref = point.get('elementRef', '')
                    if '_anchor' in element_ref:
                        # This is an anchor interaction, skip it
                        skip_interaction = True
                        break
                if skip_interaction:
                    continue

            # Get interaction type from arrowHead attribute in Point elements
            interaction_type = interaction.get('type', None)

            # If no type attribute, look for arrowHead in waypoints
            if not interaction_type:
                for point in interaction.findall('.//gpml:Point', ns):
                    arrow_head = point.get('arrowHead', None)
                    if arrow_head:
                        interaction_type = arrow_head
                        break

            # Default to Unknown if still no type found
            if not interaction_type:
                interaction_type = 'Unknown'

            self.interaction_types[interaction_type] += 1
            file_stats['interaction_types'][interaction_type] += 1

            # Normalize interaction ID (remove pathway-specific prefixes)
            normalized_id = self._normalize_interaction_id(element_id)
            if normalized_id:
                self.unique_interactions.add((interaction_type, normalized_id))

            # Get Xref
            xref = interaction.find('.//gpml:Xref', ns)
            if xref is not None:
                data_source = xref.get('dataSource', 'Unknown')
                if data_source and data_source != 'Unknown':
                    self.interaction_type_xrefs[interaction_type][data_source] += 1

    def _normalize_interaction_id(self, interaction_id):
        """Normalize interaction ID by removing pathway-specific prefixes."""
        if not interaction_id:
            return None


        import re

        # Pattern 1: Gene-protein interactions "gene_protein_X_Y_Z"
        # Keep these unique - don't normalize
        if 'gene_protein' in interaction_id:
            return interaction_id  # Keep as-is to maintain uniqueness

        # Pattern 2: "_EC-X.X.X.X-RXN_..." -> "EC-X.X.X.X-RXN"
        match = re.search(r'_?((?:EC-)?\d+\.\d+\.\d+\.\d+-RXN)', interaction_id)
        if match:
            return match.group(1)

        # Pattern 3: "..._RXN-ID_..." -> "RXN-ID"
        match = re.search(r'([A-Z0-9]+-RXN(?:-[A-Z0-9]+)?)', interaction_id)
        if match:
            return match.group(1)

        # Default: return as-is
        return interaction_id

    def _process_citations(self, root, ns):
        """Process Citation elements."""
        for citation in root.findall('.//gpml:Citation', ns):
            self.total_citations += 1

            # Get citation ID
            citation_id = citation.get('elementId', '')

            # Get Xref
            xref = citation.find('.//gpml:Xref', ns)
            if xref is not None:
                identifier = xref.get('identifier', '')
                data_source = xref.get('dataSource', '').lower()

                if citation_id:
                    self.unique_citations.add(citation_id)

                # Check for PubMed
                if data_source == 'pubmed':
                    self.citations_with_pubmed += 1
                    if citation_id:
                        self.unique_citations_with_pubmed.add(citation_id)

                # Check for DOI
                elif data_source == 'doi':
                    self.citations_with_doi += 1
                    if citation_id:
                        self.unique_citations_with_doi.add(citation_id)

                # Neither PubMed nor DOI
                else:
                    self.citations_without_either += 1
                    if citation_id:
                        self.unique_citations_without_either.add(citation_id)

    def _process_properties(self, root, ns):
        """Process Property elements from DataNodes and Interactions."""
        # Process DataNode properties
        for datanode in root.findall('.//gpml:DataNode', ns):
            for prop in datanode.findall('.//gpml:Property', ns):
                prop_key = prop.get('key', '')
                if prop_key:
                    self.property_keys[prop_key] += 1
                    self.datanode_property_keys[prop_key] += 1
                    self.total_properties += 1

        # Process Interaction properties
        for interaction in root.findall('.//gpml:Interaction', ns):
            for prop in interaction.findall('.//gpml:Property', ns):
                prop_key = prop.get('key', '')
                if prop_key:
                    self.property_keys[prop_key] += 1
                    self.interaction_property_keys[prop_key] += 1
                    self.total_properties += 1

    def analyze_directory(self, directory):
        """Analyze all GPML files in a directory."""
        gpml_files = []

        # Find all GPML files
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if filename.endswith('.gpml'):
                    gpml_files.append(os.path.join(root, filename))

        print(f"\nFound {len(gpml_files)} GPML files")
        print("Processing files...")

        # Process each file
        for i, filepath in enumerate(gpml_files, 1):
            if i % 100 == 0:
                print(f"  Progress: {i}/{len(gpml_files)} ({i*100//len(gpml_files)}%)")

            self.parse_gpml_file(filepath)

        print(f"\n+ Processed {self.files_processed} files successfully")
        if self.files_failed > 0:
            print(f"X Failed to process {self.files_failed} files")

    def print_report(self):
        """Print statistics report."""
        print("\n" + "="*80)
        print("GPML STATISTICS REPORT")
        print("="*80)
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Files analyzed: {self.files_processed}")
        print("="*80)

        # Interaction Statistics
        print("\n" + "="*80)
        print("INTERACTION STATISTICS")
        print("="*80)

        print("\nInteraction Types:")
        print(f"{'Type':<30} {'Total':<15} {'Unique':<15}")
        print("-"*60)

        # Count unique by type
        unique_by_type = defaultdict(set)
        for interaction_type, element_id in self.unique_interactions:
            unique_by_type[interaction_type].add(element_id)

        total_interactions = sum(self.interaction_types.values())
        for interaction_type in sorted(self.interaction_types.keys()):
            total = self.interaction_types[interaction_type]
            unique = len(unique_by_type[interaction_type])
            print(f"{interaction_type:<30} {total:<15,} {unique:<15,}")

        print("-"*60)
        print(f"{'TOTAL':<30} {total_interactions:<15,} {len(self.unique_interactions):<15,}")

        print("\nInteraction Xref Coverage by Type:")
        print(f"{'Type':<30} {'Data Source':<20} {'Count':<15}")
        print("-"*65)

        for interaction_type in sorted(self.interaction_type_xrefs.keys()):
            xrefs = self.interaction_type_xrefs[interaction_type]
            for data_source in sorted(xrefs.keys(), key=lambda x: xrefs[x], reverse=True):
                count = xrefs[data_source]
                print(f"{interaction_type:<30} {data_source:<20} {count:<15,}")

        if not self.interaction_type_xrefs:
            print("  (No Xrefs found in interactions)")

        # DataNode Statistics
        print("\n" + "="*80)
        print("DATANODE STATISTICS")
        print("="*80)

        print("\nDataNode Types:")
        print(f"{'Type':<30} {'Total':<15} {'Unique':<15}")
        print("-"*60)

        # Count unique by type
        unique_dn_by_type = defaultdict(set)
        for node_type, unique_id in self.unique_datanodes:
            unique_dn_by_type[node_type].add(unique_id)

        total_datanodes = sum(self.datanode_types.values())
        for node_type in sorted(self.datanode_types.keys()):
            total = self.datanode_types[node_type]
            unique = len(unique_dn_by_type[node_type])
            print(f"{node_type:<30} {total:<15,} {unique:<15,}")

        # Add protein complexes
        print(f"{'Protein Complex (Group)':<30} {self.total_complexes:<15,} {len(self.unique_complexes):<15,}")

        print("-"*60)
        total_with_complexes = total_datanodes + self.total_complexes
        unique_with_complexes = len(self.unique_datanodes) + len(self.unique_complexes)
        print(f"{'TOTAL (incl. complexes)':<30} {total_with_complexes:<15,} {unique_with_complexes:<15,}")

        print("\nDataNode Xref Coverage by Type (Unique Entities):")
        print(f"{'Type':<30} {'Data Source':<20} {'Unique':<15}")
        print("-"*65)

        for node_type in sorted(self.datanode_type_xrefs_unique.keys()):
            xrefs_unique = self.datanode_type_xrefs_unique[node_type]
            for data_source in sorted(xrefs_unique.keys(), key=lambda x: len(xrefs_unique[x]), reverse=True):
                unique_count = len(xrefs_unique[data_source])
                print(f"{node_type:<30} {data_source:<20} {unique_count:<15,}")

        # Placeholder Statistics
        print("\n" + "="*80)
        print("PLACEHOLDER STATISTICS")
        print("="*80)

        print(f"\nTotal placeholders: {self.total_placeholders:,}")
        print(f"Unique placeholders: {len(self.unique_placeholders):,}")

        if self.unique_placeholders:
            print("\nPlaceholders by Type:")
            placeholder_by_type = defaultdict(int)
            for node_type, unique_id in self.unique_placeholders:
                placeholder_by_type[node_type] += 1

            print(f"{'Type':<30} {'Count':<15}")
            print("-"*45)
            for node_type in sorted(placeholder_by_type.keys()):
                count = placeholder_by_type[node_type]
                print(f"{node_type:<30} {count:<15,}")

        # Publication Statistics
        print("\n" + "="*80)
        print("="*80)

        print(f"\nTotal citations: {self.total_citations:,}")
        print(f"Unique citations: {len(self.unique_citations):,}")

        print("\nCitation Quality:")
        print(f"{'Category':<40} {'Total':<15} {'Unique':<15} {'%':<10}")
        print("-"*70)

        pubmed_pct = (self.citations_with_pubmed / self.total_citations * 100) if self.total_citations > 0 else 0
        doi_pct = (self.citations_with_doi / self.total_citations * 100) if self.total_citations > 0 else 0
        neither_pct = (self.citations_without_either / self.total_citations * 100) if self.total_citations > 0 else 0

        print(f"{'With PubMed ID':<40} {self.citations_with_pubmed:<15,} {len(self.unique_citations_with_pubmed):<15,} {pubmed_pct:>6.1f}%")
        print(f"{'With DOI':<40} {self.citations_with_doi:<15,} {len(self.unique_citations_with_doi):<15,} {doi_pct:>6.1f}%")
        print(f"{'Without PubMed/DOI (BioCyc only)':<40} {self.citations_without_either:<15,} {len(self.unique_citations_without_either):<15,} {neither_pct:>6.1f}%")

        # Taxonomy Statistics
        print("\n" + "="*80)
        print("TAXONOMY ANNOTATION STATISTICS")
        print("="*80)

        # Pathway-level organism
        print(f"\nPathway-level organism values ({len(self.pathway_organisms)} unique):")
        for org, count in sorted(self.pathway_organism_counts.items(), key=lambda x: -x[1]):
            pct = count / self.files_processed * 100 if self.files_processed else 0
            print(f"  {org:<35} {count:>6,} files  ({pct:.1f}%)")

        # Per-type taxonomy annotation coverage
        print("\nDataNode taxonomy annotation coverage:")
        print(f"{'Type':<20} {'Total':>10} {'Annotated':>12} {'Coverage':>10}")
        print("-"*55)
        for node_type in sorted(self.total_by_type.keys()):
            total = self.total_by_type[node_type]
            annotated = self.annotated_by_type.get(node_type, 0)
            pct = annotated / total * 100 if total else 0
            print(f"{node_type:<20} {total:>10,} {annotated:>12,} {pct:>9.1f}%")

        # Top species per entity type
        for node_type in sorted(self.taxonomy_by_type.keys()):
            species_counts = self.taxonomy_by_type[node_type]
            if not species_counts:
                continue
            print(f"\nTop species — {node_type}:")
            print(f"  {'Species':<35} {'Count':>8}")
            print("  " + "-"*45)
            top = sorted(species_counts.items(), key=lambda x: -x[1])[:15]
            for species, count in top:
                print(f"  {species:<35} {count:>8,}")
            if len(species_counts) > 15:
                print(f"  ... and {len(species_counts) - 15} more species")

        # Species overlap across entity types
        species_by_type = {t: set(c.keys()) for t, c in self.taxonomy_by_type.items()}
        all_species = set().union(*species_by_type.values()) if species_by_type else set()
        overlap = {sp for sp in all_species
                   if sum(1 for s in species_by_type.values() if sp in s) > 1}
        print(f"\nUnique species across all DataNodes: {len(all_species)}")
        if overlap:
            print(f"Species appearing on multiple entity types ({len(overlap)}):")
            for sp in sorted(overlap):
                types_str = ", ".join(
                    f"{t} ({self.taxonomy_by_type[t][sp]:,})"
                    for t in sorted(species_by_type)
                    if sp in species_by_type[t]
                )
                print(f"  {sp}: {types_str}")

        # Pathway vs Reaction Split
        print("\n" + "="*80)
        print("PATHWAY vs REACTION FILE STATISTICS")
        print("="*80)

        print(f"\nPathway files: {self.pathway_files:,}")
        print(f"Reaction files: {self.reaction_files:,}")

        # Pathway-specific stats
        self._print_file_type_stats("PATHWAY FILES", self.pathway_stats)

        # Reaction-specific stats
        self._print_file_type_stats("REACTION FILES", self.reaction_stats)

        # Property Statistics
        print("\n" + "="*80)
        print("PROPERTY STATISTICS")
        print("="*80)

        print(f"\nTotal properties: {self.total_properties:,}")
        print(f"Unique property keys: {len(self.property_keys)}")

        print("\nTop 20 Property Keys (All Elements):")
        print(f"{'Property Key':<40} {'Count':<15} {'%':<10}")
        print("-"*65)
        sorted_props = sorted(self.property_keys.items(), key=lambda x: x[1], reverse=True)[:20]
        for prop_key, count in sorted_props:
            pct = (count / self.total_properties * 100) if self.total_properties > 0 else 0
            print(f"{prop_key:<40} {count:<15,} {pct:>6.1f}%")

        print("\nTop 15 DataNode Property Keys:")
        print(f"{'Property Key':<40} {'Count':<15}")
        print("-"*55)
        sorted_dn_props = sorted(self.datanode_property_keys.items(), key=lambda x: x[1], reverse=True)[:15]
        for prop_key, count in sorted_dn_props:
            print(f"{prop_key:<40} {count:<15,}")

        if self.interaction_property_keys:
            print("\nTop 10 Interaction Property Keys:")
            print(f"{'Property Key':<40} {'Count':<15}")
            print("-"*55)
            sorted_int_props = sorted(self.interaction_property_keys.items(), key=lambda x: x[1], reverse=True)[:10]
            for prop_key, count in sorted_int_props:
                print(f"{prop_key:<40} {count:<15,}")


    def _print_file_type_stats(self, title, stats):
        """Print statistics for a specific file type (pathways or reactions)."""
        print(f"\n--- {title} ---")

        # DataNode types
        total_datanodes = sum(stats['datanode_types'].values())
        unique_datanodes = len(stats['unique_datanodes'])

        print(f"\nDataNodes: {total_datanodes:,} total, {unique_datanodes:,} unique")

        # DataNode Xref coverage
        if stats['datanode_type_xrefs']:
            print("\nDataNode Xref Coverage:")
            print(f"{'Type':<20} {'Data Source':<25} {'Count':<15}")
            print("-"*60)
            for node_type in sorted(stats['datanode_type_xrefs'].keys()):
                xrefs = stats['datanode_type_xrefs'][node_type]
                for data_source in sorted(xrefs.keys(), key=lambda x: xrefs[x], reverse=True)[:5]:  # Top 5 per type
                    count = xrefs[data_source]
                    print(f"{node_type:<20} {data_source:<25} {count:<15,}")


    def save_report(self, output_file):
        """Save report to file."""
        import io
        import sys

        # Redirect stdout to string buffer
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO()

        # Print report to buffer
        self.print_report()

        # Get report content
        report_content = buffer.getvalue()

        # Restore stdout
        sys.stdout = old_stdout

        # Write to file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(report_content)

        print(f"\nReport saved: {output_file}")


def main():
    """Main execution function."""
    # Check command line arguments
    if len(sys.argv) != 2:
        print(__doc__)
        print("\nError: Missing required argument!")
        print("Usage: python analyze_gpml_stats.py <output_dir>")
        sys.exit(1)

    output_dir = sys.argv[1]

    # Validate directory
    if not os.path.exists(output_dir):
        print(f"Error: Directory does not exist: {output_dir}")
        sys.exit(1)

    print("="*80)
    print("GPML Statistics Analyzer")
    print("="*80)
    print(f"Analyzing directory: {output_dir}")

    # Create analyzer
    analyzer = GPMLAnalyzer()

    # Analyze directory
    analyzer.analyze_directory(output_dir)

    # Print report
    analyzer.print_report()

    # Save report
    report_file = os.path.join(output_dir, "GPML_STATISTICS_REPORT.txt")
    analyzer.save_report(report_file)


if __name__ == "__main__":
    main()
