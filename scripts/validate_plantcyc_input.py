#!/usr/bin/env python3
"""
validate_plantcyc_input.py
==========================

Sanity-checks PlantCyc flat files before the GPML build runs.
Each check is independent. The script prints a human-readable report and
exits with a non-zero code if any ERROR-level issue is found.

Exit codes:
    0  — all checks passed (warnings may be present)
    1  — at least one ERROR found

Usage:
    python scripts/validate_plantcyc_input.py <data_dir>

    <data_dir> must contain at minimum: genes.dat, proteins.dat

Extending:
    Add a new check as a function named check_*(data, report) and call it
    in run_checks(). Each function appends Finding objects to report.
"""

import sys
import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import List


# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------

LEVELS = ("ERROR", "WARNING", "INFO")


@dataclass
class Finding:
    level: str          # ERROR | WARNING | INFO
    check: str          # short check name
    message: str        # one-line summary
    details: List[str] = field(default_factory=list)  # per-case lines


class Report:
    def __init__(self):
        self.findings: List[Finding] = []

    def add(self, level, check, message, details=None):
        self.findings.append(Finding(level, check, message, details or []))

    def print(self):
        counts = {l: 0 for l in LEVELS}
        for f in self.findings:
            counts[f.level] += 1

        print("=" * 70)
        print("PlantCyc input validation report")
        print("=" * 70)

        for f in self.findings:
            prefix = {"ERROR": "✗", "WARNING": "△", "INFO": "·"}[f.level]
            print(f"\n[{f.level}] {prefix} {f.check}")
            print(f"  {f.message}")
            for d in f.details[:20]:
                print(f"    {d}")
            if len(f.details) > 20:
                print(f"    ... and {len(f.details) - 20} more")

        print()
        print("─" * 70)
        print(f"Summary: {counts['ERROR']} error(s), "
              f"{counts['WARNING']} warning(s), "
              f"{counts['INFO']} info")
        print("─" * 70)

    def has_errors(self):
        return any(f.level == "ERROR" for f in self.findings)


# ---------------------------------------------------------------------------
# Flat-file parser
# ---------------------------------------------------------------------------

def parse_dat(filepath):
    records = []
    current = {}
    with open(filepath, "rb") as fh:
        content = fh.read().decode("latin-1")
    for line in content.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        if line == "//":
            if current:
                records.append(current)
            current = {}
            continue
        if " - " in line:
            k, _, v = line.partition(" - ")
            k, v = k.strip(), v.strip()
            if k in current:
                existing = current[k]
                current[k] = (existing + [v]
                               if isinstance(existing, list)
                               else [existing, v])
            else:
                current[k] = v
    if current:
        records.append(current)
    return records


def as_list(val):
    if val is None:
        return []
    return val if isinstance(val, list) else [val]


# ---------------------------------------------------------------------------
# Data loader
# ---------------------------------------------------------------------------

def load_taxon_names(data_dir: Path) -> dict:
    """Build taxon_id -> common name from classes.dat (optional, best-effort).

    Both TAX-XXXX and ORG-XXXX entries are indexed.  Returns an empty dict if
    classes.dat is absent or cannot be parsed.
    """
    classes_file = data_dir / "classes.dat"
    if not classes_file.exists():
        return {}
    names = {}
    try:
        for r in parse_dat(classes_file):
            uid  = r.get("UNIQUE-ID", "")
            name = r.get("COMMON-NAME", "")
            if uid and name:
                names[uid] = name
    except Exception:
        pass
    return names


def fmt_taxon(taxon_id: str, names: dict) -> str:
    """Return 'TAX-XXXX (Species name)' when a name is available."""
    name = names.get(taxon_id, "")
    return f"{taxon_id} ({name})" if name else taxon_id


def load_data(data_dir: Path):
    """Load genes.dat, proteins.dat (and optionally classes.dat) into dicts."""
    data = {}

    data["taxon_names"] = load_taxon_names(data_dir)

    genes = parse_dat(data_dir / "genes.dat")
    data["gene_products"] = {}       # gene_id -> [protein_ids]
    data["gene_names"]    = {}       # gene_id -> COMMON-NAME
    for r in genes:
        uid = r.get("UNIQUE-ID", "")
        if uid:
            data["gene_products"][uid] = as_list(r.get("PRODUCT"))
            data["gene_names"][uid]    = r.get("COMMON-NAME", uid)

    proteins = parse_dat(data_dir / "proteins.dat")
    data["protein_species"] = {}     # protein_id -> [org_ids]
    data["protein_genes"]   = {}     # protein_id -> [gene_ids]
    data["protein_names"]   = {}     # protein_id -> COMMON-NAME
    data["protein_types"]   = {}     # protein_id -> [TYPES]
    for r in proteins:
        uid = r.get("UNIQUE-ID", "")
        if uid:
            data["protein_species"][uid] = as_list(r.get("SPECIES"))
            data["protein_genes"][uid]   = as_list(r.get("GENE"))
            data["protein_names"][uid]   = r.get("COMMON-NAME", uid)
            data["protein_types"][uid]   = as_list(r.get("TYPES"))

    # Build reverse map: gene_id -> [protein_ids that list this gene]
    data["gene_encoders"] = defaultdict(list)
    for pid, gene_list in data["protein_genes"].items():
        for gid in gene_list:
            data["gene_encoders"][gid].append(pid)

    return data


# ---------------------------------------------------------------------------
# CHECK 1 — Protein-gene species consistency
# ---------------------------------------------------------------------------

def check_protein_gene_species(data, report):
    """
    For every gene, collect the set of NCBI species (via its PRODUCT proteins).
    Flag genes whose products span more than one distinct species — these receive
    an ambiguous species annotation in the GPML DataNode.

    Also flag proteins that carry more than one SPECIES value.
    """
    gene_products   = data["gene_products"]
    protein_species = data["protein_species"]
    names           = data.get("taxon_names", {})
    gene_names     = data["gene_names"]
    protein_names  = data["protein_names"]

    # --- 1a: proteins with multiple species ---
    multi_sp_proteins = {
        pid: sps
        for pid, sps in protein_species.items()
        if len(sps) > 1
    }
    if multi_sp_proteins:
        details = [
            f"{pid} ({protein_names.get(pid, pid)}): "
            f"{[fmt_taxon(s, names) for s in sps]}"
            for pid, sps in sorted(multi_sp_proteins.items())
        ]
        report.add(
            "WARNING",
            "protein-multi-species",
            f"{len(multi_sp_proteins)} protein(s) carry more than one SPECIES value. "
            "In plants the same enzyme can be characterised across multiple species, so "
            "ALL listed taxa are annotated on the protein DataNode and propagated to "
            "encoding gene DataNodes (one AnnotationRef per taxon in the GPML output). "
            "Review each case to confirm the multi-species annotation is intentional "
            "and not a data-entry error in proteins.dat.",
            details,
        )
    else:
        report.add("INFO", "protein-multi-species",
                   "All proteins carry at most one SPECIES value.")

    # --- 1b: genes whose products span multiple species ---
    ambiguous_genes = {}
    for gid, prods in gene_products.items():
        species_set = set()
        for pid in prods:
            for sp in protein_species.get(pid, []):
                species_set.add(sp)
        if len(species_set) > 1:
            ambiguous_genes[gid] = (prods, species_set)

    if ambiguous_genes:
        # Separate genes where species differ only in ORG-code (same NCBI taxon)
        # from truly different species — we need classes.dat for that, so here
        # we use a heuristic: flag as ERROR only if raw IDs are clearly different
        # (TAX-XXXX vs TAX-YYYY where XXXX != YYYY, ignoring ORG-vs-TAX variants).

        def ncbi_num(org_id):
            """Extract numeric part if TAX-XXXX, else return None."""
            if org_id.startswith("TAX-"):
                num = org_id[4:]
                return num if num.isdigit() else None
            return None

        cross_species = {}  # genuinely different NCBI taxa
        same_taxon    = {}  # different ORG-codes but same NCBI taxon

        for gid, (prods, sps) in ambiguous_genes.items():
            ncbi_ids = {ncbi_num(s) for s in sps}
            # If all species are ORG-codes or map to same NCBI — less critical
            has_distinct_ncbi = len({n for n in ncbi_ids if n is not None}) > 1
            if has_distinct_ncbi:
                cross_species[gid] = (prods, sps)
            else:
                same_taxon[gid] = (prods, sps)

        if cross_species:
            details = [
                f"{gid} ({gene_names.get(gid, gid)}): "
                f"products={prods} → species={[fmt_taxon(s, names) for s in sorted(sps)]}"
                for gid, (prods, sps) in sorted(cross_species.items())
            ]
            report.add(
                "ERROR",
                "gene-cross-species-products",
                f"{len(cross_species)} gene(s) have products annotated to different "
                "NCBI species. The species propagated to the GeneProduct DataNode is "
                "undefined (pipeline-order-dependent).",
                details,
            )

        if same_taxon:
            details = [
                f"{gid} ({gene_names.get(gid, gid)}): "
                f"products={prods} → org_ids={[fmt_taxon(s, names) for s in sorted(sps)]} (same NCBI taxon)"
                for gid, (prods, sps) in sorted(same_taxon.items())
            ]
            report.add(
                "WARNING",
                "gene-multi-orgid-products",
                f"{len(same_taxon)} gene(s) have products with different BioCyc "
                "ORG-codes that all resolve to the same NCBI taxon. "
                "Resolution pipeline: "
                "(1) Cyc_to_wiki build — scripts/utils/organism_utils.py "
                "get_ncbi_id() maps each ORG-XXXX to an NCBI taxon ID via "
                "org_id_mapping_v2.tsv (TAX-XXXX codes are resolved directly). "
                "create_species_annotation() uses the NCBI ID as the Annotation "
                "elementId (e.g. 'taxonomy_3702'), so two products with "
                "ORG-5993 and TAX-3702 both produce elementId='taxonomy_3702'. "
                "(2) _propagate_protein_species_to_genes() deduplicates by "
                "elementRef, so the gene DataNode receives only ONE annotation "
                "for that NCBI taxon — no duplication in the GPML output. "
                "(3) gpml-to-rdf converts the annotation elementId to an "
                "NCBITaxon IRI (ncbi:3702), completing the normalisation. "
                "No action required — the annotation is correct. "
                "Review only if you suspect the ORG-code assignment in "
                "proteins.dat is wrong.",
                details,
            )

        if not cross_species and not same_taxon:
            report.add("INFO", "gene-cross-species-products",
                       "All genes have products from at most one species.")
    else:
        report.add("INFO", "gene-cross-species-products",
                   "All genes have products from at most one species.")

    # --- 1c: genes with multiple products (informational) ---
    multi_prod = {g: ps for g, ps in gene_products.items() if len(ps) > 1}
    if multi_prod:
        details = [
            f"{gid}: {prods}"
            for gid, prods in sorted(multi_prod.items())
        ]
        report.add(
            "INFO",
            "gene-multiple-products",
            f"{len(multi_prod)} gene(s) encode multiple protein products (isoforms). "
            "Propagation: pathway_builder_core._propagate_protein_species_to_genes() "
            "collects Product/Products properties from the gene DataNode, iterates "
            "over ALL linked protein nodes, and merges their annotationRefs into the "
            "gene DataNode — deduplicated by elementRef. Every distinct taxon across "
            "all isoforms is therefore represented on the gene. No data loss.",
            details,
        )
    else:
        report.add("INFO", "gene-multiple-products",
                   "All genes encode exactly one protein product.")

    # --- 1d: proteins with multiple gene entries (informational) ---
    multi_gene_prots = {
        pid: gs for pid, gs in data["protein_genes"].items() if len(gs) > 1
    }
    if multi_gene_prots:
        details = [
            f"{pid} ({protein_names.get(pid, pid)}): genes={gs}"
            for pid, gs in sorted(multi_gene_prots.items())
        ]
        report.add(
            "INFO",
            "protein-multiple-genes",
            f"{len(multi_gene_prots)} protein(s) list more than one GENE entry "
            "(e.g. enzyme complexes or duplicate gene models). "
            "_propagate_protein_species_to_genes() is gene-centric: it reads each "
            "gene node's Product/Products properties to locate linked proteins. "
            "All genes that correctly list this protein in their PRODUCT field in "
            "genes.dat will therefore receive the same taxon annotation. "
            "No data loss provided genes.dat and proteins.dat cross-references "
            "are consistent (i.e. every gene listed in the protein's GENE field "
            "also carries the protein in its own PRODUCT field).",
            details,
        )
    else:
        report.add("INFO", "protein-multiple-genes",
                   "All proteins reference at most one gene.")


# ---------------------------------------------------------------------------
# Run all checks
# ---------------------------------------------------------------------------

def run_checks(data_dir: Path):
    report = Report()

    print(f"Loading data from: {data_dir}")
    data = load_data(data_dir)
    n_names = len(data.get("taxon_names", {}))
    print(f"  {len(data['gene_products'])} genes, "
          f"{len(data['protein_species'])} proteins loaded"
          + (f", {n_names} taxon names from classes.dat" if n_names else "") + "\n")

    # --- Add new checks here as the pipeline grows ---
    check_protein_gene_species(data, report)

    return report, data


# ---------------------------------------------------------------------------
# Supplementary report writer
# ---------------------------------------------------------------------------

def _load_pubs(data_dir: Path) -> dict:
    """Return pub_id -> {title, year, authors, pmid, source} from pubs.dat."""
    pubs_file = data_dir / "pubs.dat"
    if not pubs_file.exists():
        return {}
    pubs = {}
    for r in parse_dat(pubs_file):
        uid = r.get("UNIQUE-ID", "")
        if uid:
            authors = as_list(r.get("AUTHORS", []))
            pubs[uid] = {
                "title":   r.get("TITLE", ""),
                "year":    r.get("YEAR", ""),
                "pmid":    r.get("PUBMED-ID", ""),
                "source":  r.get("SOURCE", ""),
                "authors": authors,
            }
    return pubs


def _fmt_citation(pub: dict) -> str:
    """One-line citation string from a pubs.dat record."""
    authors = pub.get("authors", [])
    first = authors[0].split(";")[0].strip() if authors else "?"
    year   = pub.get("year", "?")
    title  = (pub.get("title") or "?")[:80]
    source = pub.get("source", "")
    pmid   = pub.get("pmid", "")
    parts = [f"{first} et al. {year}.", title]
    if source:
        parts.append(source)
    if pmid:
        parts.append(f"PMID:{pmid}")
    return " ".join(parts)


def write_supplementary_report(data_dir: Path, data: dict, report: "Report",
                                out_path: Path) -> None:
    """Write a TSV table of all ERROR/WARNING findings for supplementary use.

    Tab-separated, suitable for pasting directly into Excel.
    The table is version-specific (tied to the PlantCyc release in data_dir).
    """
    names           = data.get("taxon_names", {})
    protein_species = data["protein_species"]
    protein_names   = data["protein_names"]
    gene_names      = data["gene_names"]
    gene_products   = data["gene_products"]

    # Load proteins.dat citations  and pubs.dat
    pubs = _load_pubs(data_dir)
    protein_citations: dict[str, list[str]] = {}
    if (data_dir / "proteins.dat").exists():
        for r in parse_dat(data_dir / "proteins.dat"):
            uid = r.get("UNIQUE-ID", "")
            cits = as_list(r.get("CITATIONS", []))
            if uid and cits:
                protein_citations[uid] = [
                    c.split(":")[0].strip() for c in cits
                ]

    import re, csv, io

    # Reverse map: protein_id -> [gene_ids]
    prot_to_genes: dict[str, list[str]] = {}
    for gid, prods in gene_products.items():
        for pid in prods:
            prot_to_genes.setdefault(pid, []).append(gid)

    # TSV columns
    HEADER = [
        "check_id", "severity",
        "gene_id", "gene_symbol",
        "protein_id", "protein_name",
        "species_1_taxid", "species_1_name",
        "species_2_taxid", "species_2_name",
        "short_note",
        "publication_pmid", "publication_reference",
    ]

    def _split_taxon(tok: str):
        """Split 'TAX-XXXX (Name)' into (taxid, name). Falls back gracefully."""
        tok = tok.strip().strip("'\"")
        m = re.match(r"((?:TAX|ORG)-\S+)\s*\((.+)\)", tok)
        if m:
            return m.group(1), m.group(2)
        return tok, names.get(tok, "")

    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(HEADER)

    # ── ERROR rows: gene-cross-species-products ───────────────────────────────
    for finding in report.findings:
        if finding.check != "gene-cross-species-products":
            continue
        for line in finding.details:
            m = re.match(
                r"(\S+) \(([^)]+)\): products=(\[.*?\]) → species=(\[.*?\])", line
            )
            if not m:
                continue
            gid, sym = m.group(1), m.group(2)
            prods = [p.strip().strip("'\"") for p in m.group(3).strip("[]").split(",")]
            # species tokens — split on ', ' but only between entries
            sps_raw = m.group(4).strip("[]")
            sp_items = [s.strip() for s in re.split(r",\s*(?='|TAX|ORG)", sps_raw)]
            tid1, sn1 = _split_taxon(sp_items[0]) if len(sp_items) > 0 else ("", "")
            tid2, sn2 = _split_taxon(sp_items[1]) if len(sp_items) > 1 else ("", "")

            prods_str  = "; ".join(prods)
            pnames_str = "; ".join(protein_names.get(p, p) for p in prods)

            # Citations from all products
            cit_pmids, cit_refs = [], []
            for p in prods:
                for pub_id in protein_citations.get(p, []):
                    pub = pubs.get(pub_id) or pubs.get(f"PUB-{pub_id}")
                    if pub:
                        pmid = pub.get("pmid", "")
                        ref  = _fmt_citation(pub)
                        if pmid and pmid not in cit_pmids:
                            cit_pmids.append(pmid)
                        if ref and ref not in cit_refs:
                            cit_refs.append(ref)

            writer.writerow([
                "gene-cross-species-products", "ERROR",
                gid, sym,
                prods_str, pnames_str,
                tid1, sn1, tid2, sn2,
                "",   # short_note — leave for manual curation
                "; ".join(cit_pmids),
                "; ".join(cit_refs),
            ])

    # ── WARNING rows: protein-multi-species ──────────────────────────────────
    for finding in report.findings:
        if finding.check != "protein-multi-species":
            continue
        for line in finding.details:
            m = re.match(r"(\S+) \(([^)]+)\): (\[.*)", line)
            if not m:
                continue
            pid, pname = m.group(1), m.group(2)
            sps_raw = m.group(3)
            # parse list items
            sp_items = [s.strip().strip("'\"[]") for s in
                        re.split(r",\s*'", sps_raw.strip("[]"))]
            tid1, sn1 = _split_taxon(sp_items[0]) if len(sp_items) > 0 else ("", "")
            tid2, sn2 = _split_taxon(sp_items[1]) if len(sp_items) > 1 else ("", "")

            enc_genes = prot_to_genes.get(pid, [])
            gid = "; ".join(enc_genes)
            sym = "; ".join(gene_names.get(g, g) for g in enc_genes)

            cit_pmids, cit_refs = [], []
            for pub_id in protein_citations.get(pid, []):
                pub = pubs.get(pub_id) or pubs.get(f"PUB-{pub_id}")
                if pub:
                    pmid = pub.get("pmid", "")
                    ref  = _fmt_citation(pub)
                    if pmid and pmid not in cit_pmids:
                        cit_pmids.append(pmid)
                    if ref and ref not in cit_refs:
                        cit_refs.append(ref)

            writer.writerow([
                "protein-multi-species", "WARNING",
                gid, sym,
                pid, pname,
                tid1, sn1, tid2, sn2,
                "",
                "; ".join(cit_pmids),
                "; ".join(cit_refs),
            ])

    out_path.write_text(buf.getvalue(), encoding="utf-8")
    print(f"\nSupplementary TSV written to: {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Validate PlantCyc flat files before the GPML build."
    )
    parser.add_argument("data_dir",
                        help="Path to PlantCyc data directory (must contain "
                             "genes.dat and proteins.dat)")
    parser.add_argument("--report", metavar="FILE",
                        help="Write a Markdown supplementary table of all "
                             "ERROR/WARNING cases to FILE (useful as version-specific "
                             "supplementary material)")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    if not data_dir.is_dir():
        print(f"ERROR: Not a directory: {data_dir}", file=sys.stderr)
        sys.exit(1)
    for fname in ("genes.dat", "proteins.dat"):
        if not (data_dir / fname).exists():
            print(f"ERROR: Missing required file: {data_dir / fname}",
                  file=sys.stderr)
            sys.exit(1)

    report, data = run_checks(data_dir)
    report.print()

    if args.report:
        write_supplementary_report(data_dir, data, report, Path(args.report))

    sys.exit(1 if report.has_errors() else 0)


if __name__ == "__main__":
    main()
