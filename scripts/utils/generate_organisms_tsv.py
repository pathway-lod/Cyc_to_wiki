#!/usr/bin/env python3
"""
Generate organisms.tsv from GPML files and (optionally) merge with an existing organisms.tsv.

Outputs:
  - organisms.tsv (clean TSV)
  - organisms.tsv.meta.json (provenance + counts)

Usage:
  python generate_organisms_tsv.py <gpml_dir> \
    --existing <bridgedb_organisms.tsv> \
    --output <out_tsv> \
    --meta-output <out_meta_json> \
    --bridgedb-commit <sha> \
    --bridgedb-sha256 <sha> \
    --gpml-input-dir <path> \
    --gpml-dir-name <name> \
    --generator-git-commit <sha>
"""

import xml.etree.ElementTree as ET
from pathlib import Path
import argparse
import json
from datetime import datetime, timezone


def extract_organism_info(directory: str):
    """Extract {latin_name: ncbi_taxon_id} from GPML files."""
    organisms = {}

    # Legacy fallback mapping for cases where only header organism exists
    FALLBACK_IDS = {
        'Cirsium': '41549',
        'Malus hupehensis': '106556',
        'Aveninae': '640623',
        'Rauvolfia': '4059',
        'Secale': '4549',
        'Arabidopsis thaliana': '3702',
        'Albizia': '3812'
    }

    ns = {'gpml': 'http://pathvisio.org/GPML/2021'}

    for gpml_file in Path(directory).rglob("*.gpml"):
        try:
            tree = ET.parse(gpml_file)
            root = tree.getroot()

            # 1) Pathway header organism attribute (fallback mapping only)
            org_attr = root.get('organism')
            if org_attr:
                header_orgs = [o.strip() for o in org_attr.split(',') if o.strip()]
                for org_name in header_orgs:
                    if org_name in FALLBACK_IDS:
                        organisms[org_name] = FALLBACK_IDS[org_name]

            # 2) Taxonomy annotations (preferred)
            annotations = root.findall('.//gpml:Annotation[@type="Taxonomy"]', ns)
            if not annotations:
                annotations = root.findall('.//Annotation[@type="Taxonomy"]')

            for ann in annotations:
                latin_name = ann.get('value')

                xref = ann.find('gpml:Xref', ns)
                if xref is None:
                    xref = ann.find('Xref')

                if latin_name and xref is not None:
                    identifier = xref.get('identifier')
                    datasource = xref.get('dataSource')
                    if datasource == 'NCBI Taxonomy' and identifier and identifier.isdigit():
                        organisms[latin_name] = identifier

        except Exception:
            # GPML might be malformed in rare cases; ignore and continue
            pass

    return organisms


def make_symbols_unique(organisms_list):
    """Ensure all symbols are unique in-place."""
    used = set()

    for org in organisms_list:
        sym = org.get('symbol', '')
        if not sym:
            continue

        if sym not in used:
            used.add(sym)
            continue

        genus = org.get('genus', '')
        species = org.get('species', '')

        # Try 3-letter code
        if genus and species and len(genus) >= 2:
            candidate = genus[:2].capitalize() + species[0].lower()
            if candidate not in used:
                org['symbol'] = candidate
                used.add(candidate)
                continue

        # Otherwise append numbers
        base = sym
        i = 2
        while f"{base}{i}" in used:
            i += 1
        org['symbol'] = f"{base}{i}"
        used.add(org['symbol'])


def load_existing_organisms(filepath: str):
    """Load existing organisms from TSV. Returns dict keyed by NCBI ID."""
    existing = {}
    path = Path(filepath)
    if not path.exists():
        return existing

    lines = path.read_text(encoding='utf-8').splitlines()
    if len(lines) <= 1:
        return existing

    for line in lines[1:]:
        parts = line.split('\t')
        if len(parts) < 5:
            continue
        genus, species, short_name, symbol, ncbi_id = parts[:5]
        existing[ncbi_id] = {
            'genus': genus,
            'species': species,
            'short_name': short_name,
            'symbol': symbol,
            'ncbi': ncbi_id,
            'scientific_name': short_name,
        }

    return existing


def build_new_organisms(organism_info: dict, existing_by_ncbi: dict):
    """Create new organism entries from extracted mapping, skipping existing."""
    new = []
    existing_names = {org.get('scientific_name') for org in existing_by_ncbi.values()}

    for latin_name, ncbi_id in organism_info.items():
        if ncbi_id in existing_by_ncbi:
            continue
        if latin_name in existing_names:
            continue

        parts = latin_name.split()
        genus = parts[0] if len(parts) > 0 else ''
        species = parts[1] if len(parts) > 1 else ''

        symbol = ''
        if genus and species:
            symbol = genus[0].upper() + species[0].lower()

        new.append({
            'genus': genus,
            'species': species,
            'short_name': latin_name,
            'symbol': symbol,
            'ncbi': ncbi_id,
            'scientific_name': latin_name,
        })

    return new


def write_tsv(path: str, organisms_list):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open('w', encoding='utf-8') as f:
        f.write("genus\tspecies\tshort_name\tsymbol\tncbi\n")
        for org in organisms_list:
            f.write(
                f"{org.get('genus','')}\t{org.get('species','')}\t{org.get('short_name','')}\t"
                f"{org.get('symbol','')}\t{org.get('ncbi','')}\n"
            )


def write_meta(path: str, meta: dict):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description="Generate organisms.tsv from GPML files and provenance meta.json")
    parser.add_argument('directory', help='Directory containing GPML files')
    parser.add_argument('--output', required=True, help='Output TSV file path')
    parser.add_argument('--meta-output', required=True, help='Output meta JSON file path')
    parser.add_argument('--existing', default=None, help='Existing organisms.tsv file to merge with (optional)')

    # Provenance fields (provided by wrapper script)
    parser.add_argument('--bridgedb-commit', default=None, help='BridgeDb datasources git commit SHA')
    parser.add_argument('--bridgedb-sha256', default=None, help='SHA256 of the BridgeDb organisms.tsv used')
    parser.add_argument('--gpml-input-dir', default=None, help='Full path to GPML input directory')
    parser.add_argument('--gpml-dir-name', default=None, help='Basename of the GPML directory used')
    parser.add_argument('--generator-git-commit', default=None, help='Git commit SHA of this repo when generating')
    args = parser.parse_args()

    gpml_dir = args.directory

    existing_organisms = load_existing_organisms(args.existing) if args.existing else {}
    organism_info = extract_organism_info(gpml_dir)
    new_organisms = build_new_organisms(organism_info, existing_organisms)

    all_organisms = list(existing_organisms.values()) + new_organisms
    make_symbols_unique(all_organisms)
    all_organisms.sort(key=lambda x: (x.get('genus', '').lower(), x.get('species', '').lower()))

    write_tsv(args.output, all_organisms)

    meta = {
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "generator_git_commit": args.generator_git_commit,
        "bridgedb_git_commit": args.bridgedb_commit,
        "bridgedb_file_sha256": args.bridgedb_sha256,
        "gpml": {
            "input_directory": args.gpml_input_dir,
            "input_directory_name": args.gpml_dir_name,
        },
        "counts": {
            "existing_loaded": len(existing_organisms),
            "unique_organisms_in_gpml": len(organism_info),
            "new_added": len(new_organisms),
            "total_written": len(all_organisms),
        }
    }
    write_meta(args.meta_output, meta)

    print(f"Wrote TSV : {args.output}")
    print(f"Wrote meta: {args.meta_output}")
    print(json.dumps(meta["counts"], indent=2))


if __name__ == "__main__":
    main()