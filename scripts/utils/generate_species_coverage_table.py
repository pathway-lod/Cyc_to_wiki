#!/usr/bin/env python3
"""
Generate species coverage table from PlantCyc flat files and GPML output.

Species list is derived directly from PlantCyc flat files (classes.dat,
proteins.dat, compounds.dat) — no dependency on organisms.tsv.

For each species the table reports:
  - What PlantCyc says (proteins, genes, metabolites, pathways in source data)
  - What actually ended up in the GPML files (by DataNode type and file count)

Usage:
    python scripts/utils/generate_species_coverage_table.py \
        --data-dir ../plantcyc/17.0.0/data \
        --gpml-dir output_gpml/plantcyc17.0.0-gpml2021__gitXXX__YYYYMMDD \
        --output species_coverage.tsv
"""

import argparse
import csv
import re
import os
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def strip_html(text):
    if not text:
        return text
    text = re.sub(r'<[^>]+>', '', text)
    for entity, char in [
        ('&alpha;', 'α'), ('&beta;', 'β'), ('&gamma;', 'γ'), ('&delta;', 'δ'),
        ('&epsilon;', 'ε'), ('&mu;', 'μ'), ('&omega;', 'ω'), ('&sigma;', 'σ'),
        ('&kappa;', 'κ'), ('&Lambda;', 'Λ'), ('&lambda;', 'λ'),
        ('&rarr;', '→'), ('&larr;', '←'), ('&harr;', '↔'),
        ('&amp;', '&'), ('&lt;', '<'), ('&gt;', '>'), ('&nbsp;', ' '), ('&#160;', ' '),
    ]:
        text = text.replace(entity, char)
    return text.strip()


def parse_dat_file(filepath):
    """Parse a BioCyc flat-file (.dat) into a list of record dicts.
    Handles mixed CRLF/LF line endings and latin-1 encoding."""
    records = []
    current = {}

    with open(filepath, 'rb') as fh:
        raw = fh.read().decode('latin-1')

    for line in raw.replace('\r\n', '\n').replace('\r', '\n').split('\n'):
        line = line.rstrip('\n')
        if line.startswith('#') or line.strip() == '':
            continue
        if line == '//':
            if current:
                records.append(current)
                current = {}
            continue
        if ' - ' in line:
            key, _, value = line.partition(' - ')
            key = key.strip()
            value = value.strip()
            if key in current:
                existing = current[key]
                if isinstance(existing, list):
                    existing.append(value)
                else:
                    current[key] = [existing, value]
            else:
                current[key] = value

    if current:
        records.append(current)

    return records


def _as_list(val):
    if val is None:
        return []
    return val if isinstance(val, list) else [val]


def _fmt(items, names_dict, max_show=5):
    if not items:
        return ''
    items = sorted(set(items))
    shown = [names_dict.get(i, i) for i in items[:max_show]]
    label = '; '.join(shown)
    if len(items) > max_show:
        label += f' (... {len(items)} total)'
    return label


def _fmt_ids(items, max_show=5):
    if not items:
        return ''
    items = sorted(set(items))
    shown = items[:max_show]
    label = '; '.join(shown)
    if len(items) > max_show:
        label += f' (... {len(items)} total)'
    return label


# ---------------------------------------------------------------------------
# Load organism mapping from classes.dat
# ---------------------------------------------------------------------------

def load_org_mapping(classes_dat_path):
    """
    Build org_id -> {'name': str, 'ncbi_id': str} from classes.dat.

    For ORG-XXXX: COMMON-NAME is the scientific name, ncbi_id from TYPES - TAX-XXX.
    For TAX-XXXX: COMMON-NAME is the scientific name, ncbi_id is the numeric part.
    """
    mapping = {}
    records = parse_dat_file(classes_dat_path)

    for rec in records:
        uid = rec.get('UNIQUE-ID', '').strip()
        if not uid:
            continue

        name = strip_html(rec.get('COMMON-NAME', '')) or ''

        # Derive NCBI ID
        if uid.startswith('TAX-'):
            ncbi_id = uid[4:]  # strip "TAX-"
        else:
            # For ORG-XXXX: look for TYPES - TAX-XXXX
            ncbi_id = ''
            for t in _as_list(rec.get('TYPES')):
                t = t.strip()
                if t.startswith('TAX-'):
                    ncbi_id = t[4:]
                    break

        if name:
            mapping[uid] = {'name': name, 'ncbi_id': ncbi_id}

    return mapping


# ---------------------------------------------------------------------------
# Load flat-file data
# ---------------------------------------------------------------------------

def load_proteins(proteins_file):
    """
    Returns:
        species_proteins : org_id -> [protein_ids]
        protein_names    : protein_id -> display name
        protein_records  : protein_id -> full record dict
    """
    records = parse_dat_file(proteins_file)
    species_proteins = defaultdict(list)
    protein_names = {}
    protein_records = {}

    for rec in records:
        uid = rec.get('UNIQUE-ID', '').strip()
        if not uid:
            continue
        protein_records[uid] = rec
        protein_names[uid] = strip_html(rec.get('COMMON-NAME', uid))
        for sp in _as_list(rec.get('SPECIES')):
            org_id = sp.strip()
            if org_id:
                species_proteins[org_id].append(uid)

    return species_proteins, protein_names, protein_records


def load_genes(genes_file):
    """Returns gene_product: gene_id -> [protein_ids], gene_names: gene_id -> name."""
    records = parse_dat_file(genes_file)
    gene_product = {}
    gene_names = {}
    for rec in records:
        uid = rec.get('UNIQUE-ID', '').strip()
        if not uid:
            continue
        gene_names[uid] = strip_html(rec.get('COMMON-NAME', uid))
        gene_product[uid] = _as_list(rec.get('PRODUCT'))
    return gene_product, gene_names


def load_compounds(compounds_file):
    """
    Returns:
        species_compounds : org_id -> [compound_ids]
        compound_names    : compound_id -> display name
    """
    records = parse_dat_file(compounds_file)
    species_compounds = defaultdict(list)
    compound_names = {}
    for rec in records:
        uid = rec.get('UNIQUE-ID', '').strip()
        if not uid:
            continue
        compound_names[uid] = strip_html(rec.get('COMMON-NAME', uid))
        for sp in _as_list(rec.get('SPECIES')):
            org_id = sp.strip()
            if org_id:
                species_compounds[org_id].append(uid)
    return species_compounds, compound_names


def load_enzrxns(enzrxns_file):
    """Returns enzrxn_id -> reaction_id."""
    records = parse_dat_file(enzrxns_file)
    enzrxn_reaction = {}
    for rec in records:
        uid = rec.get('UNIQUE-ID', '').strip()
        rxn = rec.get('REACTION', '')
        if uid and rxn:
            enzrxn_reaction[uid] = rxn[0] if isinstance(rxn, list) else rxn
    return enzrxn_reaction


def load_pathways(pathways_file):
    """Returns pathway_reactions: pathway_id -> [reaction_ids], pathway_names."""
    records = parse_dat_file(pathways_file)
    pathway_reactions = {}
    pathway_names = {}
    for rec in records:
        uid = rec.get('UNIQUE-ID', '').strip()
        if not uid:
            continue
        pathway_names[uid] = strip_html(rec.get('COMMON-NAME', uid))
        rxns = _as_list(rec.get('REACTION-LIST'))
        if rxns:
            pathway_reactions[uid] = rxns
    return pathway_reactions, pathway_names


# ---------------------------------------------------------------------------
# Derived mappings
# ---------------------------------------------------------------------------

def build_protein_to_reactions(protein_records, enzrxn_reaction):
    """protein_id -> set of reaction_ids (via CATALYZES → enzrxns)."""
    protein_reactions = defaultdict(set)
    for pid, rec in protein_records.items():
        for ezr in _as_list(rec.get('CATALYZES')):
            ezr = ezr.strip()
            if ezr in enzrxn_reaction:
                protein_reactions[pid].add(enzrxn_reaction[ezr])
    return protein_reactions


def build_reaction_to_pathways(pathway_reactions):
    """reaction_id -> set of pathway_ids."""
    rxn_to_pwy = defaultdict(set)
    for pwy_id, rxns in pathway_reactions.items():
        for rxn in rxns:
            rxn_to_pwy[rxn].add(pwy_id)
    return rxn_to_pwy


def build_protein_to_genes(gene_product):
    """protein_id -> set of gene_ids."""
    prot_to_gene = defaultdict(set)
    for gid, prots in gene_product.items():
        for pid in prots:
            prot_to_gene[pid].add(gid)
    return prot_to_gene


# ---------------------------------------------------------------------------
# GPML analysis
# ---------------------------------------------------------------------------

GPML_NS = 'http://pathvisio.org/GPML/2021'


def build_gpml_species_index(gpml_dir):
    """
    Scan all GPML files under gpml_dir and build a species index:

        species_value -> {
            'GeneProduct': int,   # total DataNode occurrences across all files
            'Protein':     int,
            'Metabolite':  int,
            'files':       set,   # pathway file IDs (stem of .gpml filenames)
        }

    species_value is the raw string in the Taxonomy Annotation (may be a
    scientific name or an unresolved ORG-/TAX- ID).
    """
    index = defaultdict(lambda: {
        'GeneProduct': 0, 'Protein': 0, 'Metabolite': 0, 'files': set()
    })

    gpml_files = []
    for root, _dirs, files in os.walk(gpml_dir):
        for fname in files:
            if fname.endswith('.gpml'):
                gpml_files.append(os.path.join(root, fname))

    print(f"  Scanning {len(gpml_files)} GPML files ...")

    for path in gpml_files:
        file_stem = Path(path).stem
        try:
            tree = ET.parse(path)
        except ET.ParseError:
            continue
        root_el = tree.getroot()

        # Build per-file annotation index: elementId -> {value, type}
        ann_index = {}
        for ann in root_el.iter(f'{{{GPML_NS}}}Annotation'):
            ann_id = ann.get('elementId', '')
            if ann_id:
                ann_index[ann_id] = {
                    'value': ann.get('value', ''),
                    'type':  ann.get('type', ''),
                }

        for dn in root_el.findall(f'.//{{{GPML_NS}}}DataNode'):
            node_type = dn.get('type', 'Unknown')
            if node_type not in ('GeneProduct', 'Protein', 'Metabolite'):
                continue
            for aref in dn.findall(f'{{{GPML_NS}}}AnnotationRef'):
                ref_id = aref.get('elementRef', '')
                ann = ann_index.get(ref_id, {})
                if ann.get('type', '').lower() == 'taxonomy':
                    sp_val = ann.get('value', '').strip()
                    if sp_val:
                        index[sp_val][node_type] += 1
                        index[sp_val]['files'].add(file_stem)

    return index


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate species coverage table from PlantCyc flat files and GPML output'
    )
    parser.add_argument('--data-dir', required=True,
                        help='Path to PlantCyc data directory (contains *.dat files)')
    parser.add_argument('--gpml-dir', required=True,
                        help='Path to GPML output directory (searched recursively)')
    parser.add_argument('--output', default='species_coverage.tsv',
                        help='Output TSV file path')
    args = parser.parse_args()

    data_dir = Path(args.data_dir)

    # ------------------------------------------------------------------
    # 1. Load organism name mapping from classes.dat
    # ------------------------------------------------------------------
    print("Loading classes.dat (organism name mapping) ...")
    org_mapping = load_org_mapping(data_dir / 'classes.dat')
    print(f"  {len(org_mapping)} organisms with names")

    # ------------------------------------------------------------------
    # 2. Load flat-file data
    # ------------------------------------------------------------------
    print("Loading proteins.dat ...")
    species_proteins, protein_names, protein_records = load_proteins(data_dir / 'proteins.dat')
    print(f"  {len(protein_records)} proteins across {len(species_proteins)} species")

    print("Loading genes.dat ...")
    gene_product, gene_names = load_genes(data_dir / 'genes.dat')
    print(f"  {len(gene_product)} genes")

    print("Loading compounds.dat ...")
    species_compounds, compound_names = load_compounds(data_dir / 'compounds.dat')
    print(f"  {len(compound_names)} compounds across {len(species_compounds)} species")

    print("Loading enzrxns.dat ...")
    enzrxn_reaction = load_enzrxns(data_dir / 'enzrxns.dat')

    print("Loading pathways.dat ...")
    pathway_reactions, pathway_names = load_pathways(data_dir / 'pathways.dat')

    # ------------------------------------------------------------------
    # 3. Build derived mappings
    # ------------------------------------------------------------------
    print("Building derived mappings ...")
    protein_reactions = build_protein_to_reactions(protein_records, enzrxn_reaction)
    rxn_to_pathways   = build_reaction_to_pathways(pathway_reactions)
    prot_to_genes     = build_protein_to_genes(gene_product)

    # species -> genes (via proteins)
    species_genes = defaultdict(set)
    for org_id, prot_ids in species_proteins.items():
        for pid in prot_ids:
            species_genes[org_id].update(prot_to_genes.get(pid, set()))

    # species -> pathways (via proteins → reactions → pathways)
    species_pathways = defaultdict(set)
    for org_id, prot_ids in species_proteins.items():
        for pid in prot_ids:
            for rxn in protein_reactions.get(pid, set()):
                species_pathways[org_id].update(rxn_to_pathways.get(rxn, set()))

    # ------------------------------------------------------------------
    # 4. Build full species list from flat files
    # ------------------------------------------------------------------
    all_org_ids = set(species_proteins.keys()) | set(species_compounds.keys())
    print(f"\nTotal unique species from flat files: {len(all_org_ids)}")

    # ------------------------------------------------------------------
    # 5. Scan GPML files
    # ------------------------------------------------------------------
    print(f"\nScanning GPML files in: {args.gpml_dir}")
    gpml_index = build_gpml_species_index(args.gpml_dir)
    print(f"  {len(gpml_index)} distinct species values found in GPML annotations")

    # ------------------------------------------------------------------
    # 6. Write output
    # ------------------------------------------------------------------
    fieldnames = [
        # Identity
        'org_id',
        'scientific_name',
        'ncbi_taxon_id',
        # PlantCyc source
        'in_proteins_dat',
        'n_proteins',
        'protein_ids',
        'in_compounds_dat',
        'n_compounds',
        'compound_ids',
        'n_genes',
        'gene_ids',
        'n_pathways_source',
        'pathway_ids_source',
        # GPML output
        'in_gpml',
        'gpml_as_geneproduct',
        'gpml_as_protein',
        'gpml_as_metabolite',
        'n_gpml_geneproduct_nodes',
        'n_gpml_protein_nodes',
        'n_gpml_metabolite_nodes',
        'n_gpml_pathway_files',
        'gpml_name_form',
    ]

    rows = []
    for org_id in sorted(all_org_ids):
        info = org_mapping.get(org_id, {})
        sci_name = info.get('name', org_id)
        ncbi_id  = info.get('ncbi_id', '')

        prots = sorted(set(species_proteins.get(org_id, [])))
        comps = sorted(set(species_compounds.get(org_id, [])))
        genes = sorted(species_genes.get(org_id, set()))
        pwys  = sorted(species_pathways.get(org_id, set()))

        # GPML lookup: try resolved name first, then raw org_id
        # (handles both pre-fix and post-fix GPML builds)
        gpml_by_name   = gpml_index.get(sci_name, {})
        gpml_by_raw_id = gpml_index.get(org_id, {})

        n_gp   = gpml_by_name.get('GeneProduct', 0) + gpml_by_raw_id.get('GeneProduct', 0)
        n_prot = gpml_by_name.get('Protein', 0)     + gpml_by_raw_id.get('Protein', 0)
        n_met  = gpml_by_name.get('Metabolite', 0)  + gpml_by_raw_id.get('Metabolite', 0)
        gpml_files = gpml_by_name.get('files', set()) | gpml_by_raw_id.get('files', set())

        in_gpml = n_gp + n_prot + n_met > 0

        # Determine how the species appears in GPML (resolved name / raw ID / both)
        in_by_name   = bool(gpml_by_name.get('GeneProduct', 0) + gpml_by_name.get('Protein', 0) + gpml_by_name.get('Metabolite', 0))
        in_by_raw    = bool(gpml_by_raw_id.get('GeneProduct', 0) + gpml_by_raw_id.get('Protein', 0) + gpml_by_raw_id.get('Metabolite', 0))
        if in_by_name and in_by_raw:
            gpml_name_form = 'both'
        elif in_by_name:
            gpml_name_form = 'resolved_name'
        elif in_by_raw:
            gpml_name_form = 'raw_id'
        else:
            gpml_name_form = 'absent'

        rows.append({
            'org_id':                   org_id,
            'scientific_name':          sci_name,
            'ncbi_taxon_id':            ncbi_id,
            'in_proteins_dat':          'yes' if prots else 'no',
            'n_proteins':               len(prots),
            'protein_ids':              _fmt_ids(prots),
            'in_compounds_dat':         'yes' if comps else 'no',
            'n_compounds':              len(comps),
            'compound_ids':             _fmt_ids(comps),
            'n_genes':                  len(genes),
            'gene_ids':                 _fmt_ids(genes),
            'n_pathways_source':        len(pwys),
            'pathway_ids_source':       _fmt_ids(pwys),
            'in_gpml':                  'yes' if in_gpml else 'no',
            'gpml_as_geneproduct':      'yes' if n_gp   > 0 else 'no',
            'gpml_as_protein':          'yes' if n_prot > 0 else 'no',
            'gpml_as_metabolite':       'yes' if n_met  > 0 else 'no',
            'n_gpml_geneproduct_nodes': n_gp,
            'n_gpml_protein_nodes':     n_prot,
            'n_gpml_metabolite_nodes':  n_met,
            'n_gpml_pathway_files':     len(gpml_files),
            'gpml_name_form':           gpml_name_form,
        })

    # Sort: in_gpml first, then by scientific name
    rows.sort(key=lambda r: (r['in_gpml'] == 'no', r['scientific_name']))

    print(f"\nWriting {len(rows)} rows to {args.output} ...")
    with open(args.output, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_total       = len(rows)
    n_in_gpml     = sum(1 for r in rows if r['in_gpml'] == 'yes')
    n_as_gp       = sum(1 for r in rows if r['gpml_as_geneproduct'] == 'yes')
    n_as_prot     = sum(1 for r in rows if r['gpml_as_protein'] == 'yes')
    n_as_met      = sum(1 for r in rows if r['gpml_as_metabolite'] == 'yes')
    n_resolved    = sum(1 for r in rows if r['gpml_name_form'] == 'resolved_name')
    n_raw_id      = sum(1 for r in rows if r['gpml_name_form'] == 'raw_id')
    n_both        = sum(1 for r in rows if r['gpml_name_form'] == 'both')
    n_absent      = sum(1 for r in rows if r['gpml_name_form'] == 'absent')
    n_with_prot   = sum(1 for r in rows if r['in_proteins_dat'] == 'yes')
    n_with_cpd    = sum(1 for r in rows if r['in_compounds_dat'] == 'yes')

    print(f"\n{'='*55}")
    print(f"Species coverage summary ({n_total} species from PlantCyc)")
    print(f"{'='*55}")
    print(f"  Source (PlantCyc flat files):")
    print(f"    In proteins.dat              : {n_with_prot}")
    print(f"    In compounds.dat             : {n_with_cpd}")
    print(f"  GPML output:")
    print(f"    In GPML (any node type)      : {n_in_gpml}")
    print(f"    As GeneProduct DataNode      : {n_as_gp}")
    print(f"    As Protein DataNode          : {n_as_prot}")
    print(f"    As Metabolite DataNode       : {n_as_met}")
    print(f"  Name form in GPML annotations:")
    print(f"    Resolved scientific name     : {n_resolved}")
    print(f"    Raw ORG-/TAX- ID             : {n_raw_id}")
    print(f"    Both (resolved + raw)        : {n_both}")
    print(f"    Absent from GPML             : {n_absent}")
    print(f"{'='*55}")
    print(f"Output: {args.output}")


if __name__ == '__main__':
    main()
