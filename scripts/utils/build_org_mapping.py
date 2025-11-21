1#!/usr/bin/env python3
"""
Build mapping from ORG-ID to NCBI taxonomy ID and Latin name from BioCyc dat files.
"""

import re
import argparse
from pathlib import Path


def parse_species_dat(species_dat_path):
    """
    Parse BioCyc dat file to extract organism info mapping.
    Handles both ORG-ID and TAX-ID entries.

    Returns:
        dict: Mapping of ORG-ID/TAX-ID to {'latin_name': str, 'ncbi_id': str}
    """
    mapping = {}

    with open(species_dat_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Split into records by double newline
    records = content.split('//\n')

    for record in records:
        organism_id = None
        latin_name = None
        ncbi_id = None

        for line in record.split('\n'):
            line = line.strip()

            # Extract ORG-ID or TAX-ID
            if line.startswith('UNIQUE-ID - ORG'):
                organism_id = line.split(' - ')[1].strip()
            elif line.startswith('UNIQUE-ID - TAX-'):
                organism_id = line.split(' - ')[1].strip()
                # Extract the numeric part for NCBI ID
                tax_num = organism_id.replace('TAX-', '')
                if tax_num.isdigit():
                    ncbi_id = tax_num

            # Extract Latin/scientific name
            elif line.startswith('COMMON-NAME - '):
                latin_name = line.split(' - ', 1)[1].strip()

            # Extract NCBI taxonomy ID from TYPES field (e.g., TYPES - TAX-3702)
            elif line.startswith('TYPES - TAX-'):
                tax_part = line.split(' - ')[1].strip()
                tax_id = tax_part.replace('TAX-', '')
                if tax_id.isdigit() and not ncbi_id:  # Only set if not already set
                    ncbi_id = tax_id

            # Also extract NCBI taxonomy ID from DBLINKS (fallback)
            elif line.startswith('DBLINKS - (NCBI-TAXONOMY-DB'):
                match = re.search(r'"(\d+)"', line)
                if match and not ncbi_id:  # Only set if not already set
                    ncbi_id = match.group(1)

        # Store mapping if we have an organism ID
        if organism_id:
            mapping[organism_id] = {
                'latin_name': latin_name or '',
                'ncbi_id': ncbi_id or ''
            }

    return mapping


def main():
    parser = argparse.ArgumentParser(description='Build ORG-ID mapping from BioCyc dat files')
    parser.add_argument('data_dir', help='Path to BioCyc data directory')
    parser.add_argument('--output', default='org_mapping.tsv', help='Output TSV file')
    args = parser.parse_args()

    # Parse both classes.dat and species.dat
    mapping = {}

    classes_dat = Path(args.data_dir) / 'classes.dat'
    if classes_dat.exists():
        print(f"Parsing {classes_dat}...")
        mapping.update(parse_species_dat(classes_dat))

    species_dat = Path(args.data_dir) / 'species.dat'
    if species_dat.exists():
        print(f"Parsing {species_dat}...")
        mapping.update(parse_species_dat(species_dat))

    print(f"Found {len(mapping)} ORG-ID entries")

    # Write TSV file
    print(f"Writing mapping to {args.output}...")
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write("org_id\tlatin_name\tncbi_id\n")
        for org_id, info in sorted(mapping.items()):
            f.write(f"{org_id}\t{info['latin_name']}\t{info['ncbi_id']}\n")

    print(f"Done! Generated {args.output}")

    # Print some stats
    with_ncbi = sum(1 for info in mapping.values() if info['ncbi_id'])
    with_name = sum(1 for info in mapping.values() if info['latin_name'])
    print(f"\nStats:")
    print(f"  Total ORG-IDs: {len(mapping)}")
    print(f"  With Latin names: {with_name}")
    print(f"  With NCBI IDs: {with_ncbi}")


if __name__ == "__main__":
    main()
