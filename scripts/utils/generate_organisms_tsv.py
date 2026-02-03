#!/usr/bin/env python3
"""
Extract organism names from GPML files and generate organisms.tsv.
Uses org_id_mapping_v2.tsv to look up NCBI IDs for all organisms.
Usage: python generate_organisms_tsv.py <directory> [--output FILE] [--existing FILE] [--mapping FILE]

[--existing FILE] : the existing organisms.tsv file from BridgeDB 
[--mapping FILE] : organism_id_mapping_v2.tsv
"""

import xml.etree.ElementTree as ET
from pathlib import Path
import argparse

def load_organism_mapping(mapping_file='org_id_mapping_v2.tsv'):
    """Load organism mapping to get NCBI IDs for Latin names."""
    latin_to_ncbi = {}

    if not Path(mapping_file).exists():
        print(f"Warning: {mapping_file} not found. Will query NCBI for all organisms.")
        return latin_to_ncbi

    with open(mapping_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('org_id'):  # Skip header
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                latin_name = parts[1]
                ncbi_id = parts[2]
                if latin_name and ncbi_id:
                    latin_to_ncbi[latin_name] = ncbi_id

    print(f"Loaded {len(latin_to_ncbi)} organism mappings from {mapping_file}")
    return latin_to_ncbi

def extract_organism_names(directory):
    """Extract all unique organism Latin names from GPML files."""
    organisms = set()
    for gpml_file in Path(directory).rglob("*.gpml"):
        try:
            root = ET.parse(gpml_file).getroot()
            # Get organism attribute
            organism_text = root.get('organism', '')

            # Also check Property elements
            for prop in root.findall('.//{http://pathvisio.org/GPML/2021}Property'):
                key = prop.get('key', '')
                if key == 'Organism' or key.startswith('TaxonomicRange'):
                    organism_text += ', ' + prop.get('value', '')

            # Extract organism names (comma-separated)
            for part in organism_text.split(','):
                name = part.strip()
                # Skip empty, skip TAX-IDs, skip ORG-IDs
                if name and not name.startswith('TAX-') and not name.startswith('ORG-'):
                    organisms.add(name)
        except:
            pass

    return sorted(organisms)



def make_symbols_unique(organisms):
    """Ensure all symbols are unique."""
    used_symbols = set()
    for org in organisms:
        # Skip empty symbols
        if not org['symbol']:
            continue

        original = org['symbol']
        if original in used_symbols:
            # Try 3-letter code first
            if org['genus'] and org['species'] and len(org['genus']) >= 2:
                new_symbol = org['genus'][:2].capitalize() + org['species'][0].lower()
                if new_symbol not in used_symbols:
                    org['symbol'] = new_symbol
                    used_symbols.add(new_symbol)
                    continue

            # Otherwise add numbers
            counter = 2
            while f"{original}{counter}" in used_symbols:
                counter += 1
            org['symbol'] = f"{original}{counter}"
            used_symbols.add(org['symbol'])
        else:
            used_symbols.add(original)

def load_existing_organisms(filepath):
    """Load existing organisms from TSV file."""
    organisms = {}
    if not Path(filepath).exists():
        return organisms

    print(f"Loading existing organisms from {filepath}...")
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        if len(lines) > 1:  # Skip header
            for line in lines[1:]:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    ncbi_id = parts[4]
                    organisms[ncbi_id] = {
                        'genus': parts[0],
                        'species': parts[1],
                        'short_name': parts[2],
                        'symbol': parts[3],
                        'ncbi': ncbi_id,
                        'scientific_name': parts[2]  # Use short_name as fallback
                    }
    print(f"  Loaded {len(organisms)} existing organisms")
    return organisms

def main():
    parser = argparse.ArgumentParser(
        description='Generate organisms.tsv from GPML files with Latin names',
        usage='%(prog)s <directory> [OPTIONS]'
    )
    parser.add_argument('directory',
                        help='Directory containing GPML files')
    parser.add_argument('--output', default='organisms.tsv', help='Output TSV file (default: organisms.tsv)')
    parser.add_argument('--existing', default=None, help='Existing organisms.tsv file to merge with (optional)')
    parser.add_argument('--mapping', default='org_id_mapping_v2.tsv', help='Organism mapping file (default: org_id_mapping_v2.tsv)')
    args = parser.parse_args()

    # Load organism mapping (Latin name -> NCBI ID)
    latin_to_ncbi = load_organism_mapping(args.mapping)

    # Load existing organisms if specified
    existing_organisms = {}
    if args.existing:
        existing_organisms = load_existing_organisms(args.existing)

    # Extract organism Latin names from GPML files
    print(f"\nExtracting organism names from {args.directory}...")
    organism_names = extract_organism_names(args.directory)
    print(f"Found {len(organism_names)} unique organisms in GPML files")

    # Build organism list
    organisms = []
    not_found = []

    for latin_name in organism_names:
        # Skip if we already have it
        if any(org['scientific_name'] == latin_name for org in existing_organisms.values()):
            continue

        # Check if we have the NCBI ID in our mapping
        if latin_name in latin_to_ncbi:
            ncbi_id = latin_to_ncbi[latin_name]
            name_parts = latin_name.split()
            organisms.append({
                'genus': name_parts[0] if name_parts else '',
                'species': name_parts[1] if len(name_parts) > 1 else '',
                'short_name': latin_name,
                'symbol': name_parts[0][:1].upper() + name_parts[1][:1].lower() if len(name_parts) > 1 else '',
                'ncbi': ncbi_id,
                'scientific_name': latin_name
            })
        else:
            # Not found in mapping
            not_found.append(latin_name)

    print(f"  {len(organisms)} organisms resolved from mapping")
    if not_found:
        print(f"  Warning: {len(not_found)} organisms not found in mapping:")
        for name in not_found[:10]:  # Show first 10
            print(f"    - {name}")
        if len(not_found) > 10:
            print(f"    ... and {len(not_found) - 10} more")

    # Merge with existing organisms
    all_organisms = list(existing_organisms.values()) + organisms
    make_symbols_unique(all_organisms)

    # Sort by genus, species
    all_organisms.sort(key=lambda x: (x['genus'].lower(), x['species'].lower()))

    # Write output
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write("genus\tspecies\tshort_name\tsymbol\tncbi\n")
        for org in all_organisms:
            f.write(f"{org['genus']}\t{org['species']}\t{org['short_name']}\t{org['symbol']}\t{org['ncbi']}\n")

    print(f"\nGenerated {args.output} with {len(all_organisms)} organisms")
    print(f"  Existing: {len(existing_organisms)}")
    print(f"  From mapping: {len(organisms)}")

if __name__ == "__main__":
    main()
