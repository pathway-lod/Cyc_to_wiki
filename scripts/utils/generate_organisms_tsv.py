#!/usr/bin/env python3
"""
Extract organism names and NCBI IDs from GPML files (via Annotations) and generate organisms.tsv.
Usage: python generate_organisms_tsv.py <directory> [--output FILE] [--existing FILE]
"""

import xml.etree.ElementTree as ET
from pathlib import Path
import argparse
import re

def extract_organism_info(directory):
    """
    Extract unique organism info (Latin name, NCBI ID) from GPML files.
    Looks for <Annotation type="Taxonomy"> elements and the 'organism' attribute in the header.
    
    Returns:
        dict: {latin_name: ncbi_id}
    """
    organisms = {}
    
    # Hardcoded fallbacks for organisms because i had some errors
    # Derived from classes.dat
    FALLBACK_IDS = {
        'Cirsium': '41549',
        'Malus hupehensis': '106556',
        'Aveninae': '640623',
        'Rauvolfia': '4059',
        'Secale': '4549',
        'Arabidopsis thaliana': '3702',
        'Albizia': '3812'
    }

    # XML namespaces
    ns = {'gpml': 'http://pathvisio.org/GPML/2021'}
    
    print(f"Scanning GPML files in {directory}...")
    
    for gpml_file in Path(directory).rglob("*.gpml"):
        try:
            tree = ET.parse(gpml_file)
            root = tree.getroot()
            
            # 1. Check 'organism' attribute in Pathway tag (header)
            org_attr = root.get('organism')
            if org_attr:
                # Handle comma separated list
                header_orgs = [o.strip() for o in org_attr.split(',') if o.strip()]
                for org_name in header_orgs:
                    # Try to resolve ID from fallback
                    if org_name in FALLBACK_IDS:
                        organisms[org_name] = FALLBACK_IDS[org_name]

            # 2. Find all Taxonomy annotations
            annotations = root.findall('.//gpml:Annotation[@type="Taxonomy"]', ns)
            if not annotations:
                annotations = root.findall('.//Annotation[@type="Taxonomy"]')
                
            for ann in annotations:
                latin_name = ann.get('value')
                
                # Find Xref child
                xref = ann.find('gpml:Xref', ns)
                if xref is None:
                    xref = ann.find('Xref')
                    
                if latin_name and xref is not None:
                    identifier = xref.get('identifier')
                    datasource = xref.get('dataSource')
                    
                    # Check if it's a valid NCBI Taxon ID (numeric)
                    if datasource == 'NCBI Taxonomy' and identifier and identifier.isdigit():
                        organisms[latin_name] = identifier
                        
        except Exception as e:
            # print(f"Error parsing {gpml_file}: {e}")
            pass

    return organisms

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
    # Removed mapping argument
    args = parser.parse_args()

    # Load existing organisms if specified
    existing_organisms = {}
    if args.existing:
        existing_organisms = load_existing_organisms(args.existing)

    # Extract organism info from GPML files
    organism_info = extract_organism_info(args.directory)
    print(f"Found {len(organism_info)} unique organisms with NCBI IDs in GPML files")

    # Build organism list
    new_organisms = []
    
    for latin_name, ncbi_id in organism_info.items():
        # Skip if we already have this NCBI ID
        if ncbi_id in existing_organisms:
            continue
            
        # Also check if we have the name but under a different ID (unlikely but safe)
        if any(org['scientific_name'] == latin_name for org in existing_organisms.values()):
            continue

        name_parts = latin_name.split()
        new_organisms.append({
            'genus': name_parts[0] if name_parts else '',
            'species': name_parts[1] if len(name_parts) > 1 else '',
            'short_name': latin_name,
            'symbol': name_parts[0][:1].upper() + name_parts[1][:1].lower() if len(name_parts) > 1 else '',
            'ncbi': ncbi_id,
            'scientific_name': latin_name
        })

    print(f"  Adding {len(new_organisms)} new organisms")

    # Merge with existing organisms
    all_organisms = list(existing_organisms.values()) + new_organisms
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
    print(f"  New: {len(new_organisms)}")

if __name__ == "__main__":
    main()
