# Release Notes

This file documents each data release and provides step-by-step instructions for preparing a new one.

---

## How to prepare a release

### 1. Finish and validate the build

Run the pipeline on the target PlantCyc version:

```bash
conda activate cyc_2_wiki
./scripts/run_pipeline.sh
```

Validate all generated GPML files — there should be **zero errors** before releasing:

```bash
python scripts/test_gpml_files.py output_gpml/<run-directory>
```

### 2. Clean up old output runs

Keep only the single output directory you intend to release. Delete all intermediate or test runs from `output_gpml/` so the repository only contains the release data.

### 3. Merge to main

Open a pull request from your feature branch into `main` and merge once the build is clean.

### 4. Tag the release

From `main`, after merging:

```bash
git tag -a plantcycX.Y.Z-gpml2021-vN -m "PlantCyc X.Y.Z to GPML 2021 (vN)"
git push origin plantcycX.Y.Z-gpml2021-vN
```

The version suffix (`v1`, `v2`, …) is set in `scripts/config.env` (`RELEASE_VERSION`). Increment it whenever you re-release the same PlantCyc+GPML combination with updated code. The pipeline script uses it automatically in both the tag name and the output directory name.

The tag triggers automatic archiving on Zenodo (configured via `.zenodo.json` / GitHub integration).

### 5. Draft the GitHub release

Go to the repository's **Releases** page on GitHub, select the new tag, and fill in:
- **Title**: `PlantCyc X.Y.Z → GPML2021`
- **Description**: copy the relevant section from this file

### 6. Update this file

Add a new release section below (copy the template from the previous release) and update:
- **Code commit** — the git short SHA of the final release build
- **Build date** — today's date
- **Zenodo DOI** — once archived
- **What's new** — bullet points summarising code changes since the last release
- **Taxonomy annotation coverage** — numbers from `GPML_STATISTICS_REPORT.txt` in the release output directory

---

## Releases

---

### plantcyc17.0.0-gpml2021

**PlantCyc version:** 17.0.0  
**GPML snapshot:** GPML2021  
**Code commit:** `e2a72d55`  (`feat/species-plants` → `main`)  
**Build date:** 2026-05-15  
**Zenodo DOI:** [10.5281/zenodo.18404067](https://doi.org/10.5281/zenodo.18404067)

#### What's new

- **Viridiplantae as pathway-level organism** — All GPML pathways now carry `organism="Viridiplantae"` (NCBI Taxonomy 33090) with a structured `<Annotation type="Taxonomy">` element, replacing the previous multi-species comma-separated string that caused XSD validation errors.
- **Per-entity taxonomy annotations** — Proteins, compounds, and genes with species information now carry individual `<AnnotationRef>` → `<Annotation type="Taxonomy">` elements with NCBI Taxonomy IDs. Genes inherit species from their linked protein product (`PRODUCT` field) since `genes.dat` has no `SPECIES` field.
- **Fixed `MISSING_GROUP_REF` validation errors** — Protein monomer DataNodes now correctly reference their parent complex group using the sanitized element ID.
- **Fixed XSD: pathway-level `<AnnotationRef>` not allowed** — The GPML2021 XSD only permits `<CitationRef>` and `<EvidenceRef>` as direct children of `<Pathway>`; the Viridiplantae reference is captured via `organism=` attribute and `<Annotation>` element instead.
- **Zero validation errors** — All 2,478 GPML files pass both local and XSD schema validation.

#### Build statistics

| Category | Count |
|---|---|
| Pathway files | 1,162 |
| Reaction files | 1,316 |
| **Total GPML files** | **2,478** |
| Valid files (0 errors) | 2,478 (100%) |
| Validation warnings | 1,145 |

**DataNodes**

| Type | Total | Unique |
|---|---|---|
| Metabolite | 23,449 | 4,843 |
| Protein | 10,715 | 3,227 |
| GeneProduct | 8,259 | 2,620 |
| Protein Complex (Group) | 1,237 | 400 |

**Interactions**

| Type | Total | Unique |
|---|---|---|
| TranscriptionTranslation | 8,280 | 132 |
| Catalysis | 13,853 | 7,086 |
| Conversion | 7,772 | 5,490 |
| Inhibition | 3,486 | 1,467 |
| Stimulation | 865 | 383 |

**Taxonomy annotation coverage**

| Entity type | Total | Annotated | Coverage |
|---|---|---|---|
| GeneProduct | 8,259 | 8,148 | 98.7% |
| Protein | 10,715 | 10,583 | 98.8% |
| Metabolite | 23,449 | 33 | 0.1% |

Unique species across all DataNodes: **422**  
Top species: *Arabidopsis thaliana* (4,117), *Glycine max* (587), *Lotus japonicus* (463), *Zea mays* (336), *Chlamydomonas reinhardtii* (318)
