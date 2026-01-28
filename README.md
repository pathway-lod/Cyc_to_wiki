# PlantCyc to Wikidata repository 

This repository contains the pipeline to transform data from Biopax format (common in PlantCyc and Biocyc repositories) to GPML2021 format [(Graphical Pathway Markup Language)](https://pathvisio.org/documentation/GPML). 

## Installation

```bash
git clone git@github.com:pathway-lod/Cyc_to_wiki.git
cd Cyc_to_wiki
```

## Usage

### Basic Command

```bash
python build_pathways.py <data_dir> <output_dir> [options]
```

### Options

- `--include-reactions`: Also build single reaction files for unused reactions
- `--layout grid`: Use grid layout (default) - fast
- `--layout forceatlas2`: Use ForceAtlas2 force-directed layout

**Examples**:
```bash
# Build pathways only with default grid layout
python build_pathways.py ./data ./output

# Build pathways AND single reactions
python build_pathways.py ./data ./output --include-reactions

# Build with ForceAtlas2 layout 
python build_pathways.py ./data ./output --layout forceatlas2

# Combine options
python build_pathways.py ./data ./output --include-reactions --layout forceatlas2
```

### Input Files

Place all PlantCyc `.dat` files in a directory.

### Output

- **Individual pathway files**: `[output_dir]/biocyc_pathways_[timestamp]/individual_pathways/*.gpml`
- **Single reaction files** (if `--include-reactions`): `individual_reactions/*.gpml`
- **Analysis report**: `GPML_STATISTICS_REPORT.txt` - Generated automatically by the analysis script 

## Project Structure

```
project/
├── scripts/
│   ├── parsing_functions/
│   │   └── parsing_utils.py           # Parse BioCyc flatfile format
│   ├── data_structure/
│   │   └── wiki_data_structure.py     # GPML2021 dataclasses (Pathway, DataNode, Interaction, etc.)
│   ├── build_functions/
│   │   ├── pathway_builder_core.py    # Main pathway building logic (includes layout)
│   │   ├── general_pathwaybuilder.py  # Entry point script
│   │   ├── build_*_data_nodes.py      # Entity builders (genes, proteins, compounds)
│   │   ├── build_*_interaction.py     # Interaction builders (reactions, regulations)
│   │   └── citation_manager.py        # Publication and citation management
│   └── object2gmpl/
│       └── gpml_writer.py             # Convert Python objects to GPML XML
└── build_pathways.py                   # Main entry point
```

---

## How It Works: Pathway Building Methodology

The converter follows a bottom-up data collection approach, starting from a pathway ID and recursively gathering all related components.

### High-Level Process

```
Pathway ID
    ↓
1. Retrieve Pathway Record & Expand Sub-Pathways
    ↓
2. Collect Reactions from Expanded List (if `--include-reactions`, individual reactions will be retrieved separately and individual_reactions GPML files will be built for them)
    ↓
3. Retrieve Compounds from Reactions
    ↓
4. Retrieve Proteins/Enzymes that Catalyze Reactions
    ↓
5. Retrieve Genes that Encode Proteins
    ↓
6. Retrieve Regulations for Reactions
    ↓
7. Collect Primary Compound Information (REACTION-LAYOUT)
    ↓
8. Calculate Positions for All Components
    ↓
9. Create DataNodes (Genes, Proteins, Compounds)
    ↓
10. Create Interactions (Gene→Protein, Protein→Reaction, Compound↔Reaction, Regulator→Reaction)
    ↓
11. Gather Citations for All Elements
    ↓
GPML Pathway (XML)
```

### Steps Explained

#### 1. Pathway Expansion and Sub-Pathways

BioCyc pathways can reference other pathways in their `REACTION-LIST` field. The converter recursively expands these references to get the complete set of reactions.

**Example**:
```
ARGDEGRAD-PWY has REACTION-LIST: [ARGININE-DEIMINASE-RXN, CITRULLINE-DEG-PWY]
                                                              ↓ (sub-pathway)
CITRULLINE-DEG-PWY has REACTION-LIST: [RXN0-5222, RXN-14196, ORNCARBAMTRANSFER-RXN]

Final expanded list: [ARGININE-DEIMINASE-RXN, RXN0-5222, RXN-14196, ORNCARBAMTRANSFER-RXN]
```

#### 2. Entity Collection (Compounds, Proteins, Genes)

From the expanded reaction list, the converter collects:
- **Compounds**: All reactants and products from reaction equations
- **Proteins/Enzymes**: All catalysts from `enzrxns.dat` (bridge between proteins and reactions)
- **Genes**: All genes that encode proteins (from `PRODUCT` field in genes.dat)


#### 3. Regulation Collection

Regulations specify how compounds or proteins **activate or inhibit** reactions. The converter:
- Parses `regulation.dat` to find regulations targeting pathway reactions
- Determines regulation mode: `+` (activation) → Stimulation arrow, `-` (inhibition) → Inhibition arrow
- If the regulation mode is not specified, the regulation is not added 
- Adds regulators to the pathway, even if they weren't part of the reaction flow

#### 4. Primary Compound Information (REACTION-LAYOUT)

BioCyc provides layout hints in the `REACTION-LAYOUT` field to indicate which compounds are the "main" reactants/products for visualization.

Examples of hints: 

**Format**:
```
(RXN-ID (:LEFT-PRIMARIES CPD1 CPD2) (:DIRECTION :L2R) (:RIGHT-PRIMARIES CPD3))
```

- `LEFT-PRIMARIES`: Compounds on the left side of the chemical equation (typically reactants)
- `RIGHT-PRIMARIES`: Compounds on the right side of the chemical equation (typically products)
- `DIRECTION`: `:L2R` (left-to-right) or `:R2L` (right-to-left) - indicates **physiological flow direction**

**Important Design Choice**: The `:DIRECTION` field controls which way the arrow points but does NOT change the underlying chemical equation.

- **L2R (left-to-right)**: Arrow points from LEFT-PRIMARIES → RIGHT-PRIMARIES (normal flow)
- **R2L (right-to-left)**: Arrow points from RIGHT-PRIMARIES → LEFT-PRIMARIES (reversed physiological flow)

**Technical Implementation**: When `DIRECTION` is `:R2L`, the converter:
1. Sets arrow start to RIGHT-PRIMARIES (products become visual start)
2. Sets arrow end to LEFT-PRIMARIES (reactants become visual end)
3. Searches for arrow start compounds in BOTH products and reactants lists (critical for R2L)
4. Searches for arrow end compounds in BOTH reactants and products lists (critical for R2L)


**Visual**:
```
Direction :L2R (normal)
FRUCTOSE-1,6-BP ────────○────────> GAP + DHAP
   (reactant)                     (products)
   Arrow flows left → right

Direction :R2L (reversed physiological flow)
FRUCTOSE-1,6-BP <──────○──────── GAP + DHAP
   (reactant)                     (products)
   Arrow flows right → left (reversed to show physiological direction)
```

**Search Strategy**: The converter searches for `REACTION-LAYOUT` across all pathway records. This handles cases where sub-pathways have more layout information than parent pathways, to ensure that layout information is kept in the parent pathway.

#### 5. Layout Calculation

The converter supports two layout algorithms, configurable via `--layout` option:

##### 5a. Grid Layout (Default)

A three-layer hierarchical grid layout organizes entities by type:

```
┌─────────────────────────────────────────┐
│  Layer 1: GENES (top)                   │
│  Grid: 3+ columns, 120px spacing        │
└─────────────────────────────────────────┘
           ↓ 80px padding
┌─────────────────────────────────────────┐
│  Layer 2: PROTEINS (middle)             │
│  Grid: 3+ columns, 150px spacing        │
└─────────────────────────────────────────┘
           ↓ 80px padding
┌─────────────────────────────────────────┐
│  Layer 3: COMPOUNDS (bottom)            │
│  Grid: 3+ columns, 150px spacing        │
└─────────────────────────────────────────┘
```

**Pros**: Fast, predictable, hierarchical view of gene→protein→compound flow
**Cons**: May have more edge crossings in complex pathways

##### 5b. ForceAtlas2 Layout

Force-directed graph layout that automatically positions nodes to minimize edge crossings and overlap. Using NetworkX module. 

Other layouts can be implemented. 


#### 6. Interaction Creation

The converter creates four types of interactions:

**A. Gene → Protein** (Transcription/Translation)
- Arrow from gene to protein with `TranscriptionTranslation` arrow head
- Connects to all monomers if protein is a complex

**B. Protein → Reaction** (Catalysis)
- Uses a **central anchor pattern**:  
```
Main Reactant ────────○──────── Main Product
                      │
                   Enzyme
```
- Main reactant/product on the main line (from `REACTION-LAYOUT`)
- Additional reactants/products connect to anchor `─○─`
- Enzyme connects to anchor from top with `Catalysis` arrow head

**C. Compound ↔ Reaction** (Conversion)
- Additional reactants: compound → anchor (no arrow)
- Additional products: anchor → compound (with `Conversion` arrow head)

**D. Regulator → Reaction** (Regulation)
- Arrow from regulator to reaction anchor
- Arrow type based on mode: `Stimulation` (+) or `Inhibition` (-)

#### 7. Citation Management

The converter extracts citations from:
- `CITATIONS` field in records
- Inline citations in `COMMENT` fields (e.g., `|CITS: [12345]|`)

**Citation Priority** (in order of preference):
1. **PubMed ID** → Creates `<Xref dataSource="pubmed" identifier="12345678" />`
   - Standard database identifier that PathVisio/WikiPathways can use to link to PubMed
2. **DOI** → Creates `<Xref dataSource="doi" identifier="12345678" />`
    - When no PubMed ID is available but DOI is, DOI is set as main Xref 
4. **BioCyc ID** → Creates `<Xref dataSource="BioCyc" identifier="UNIQUE-ID | Title (Authors, Year)" />` with full bibliographic information:
   - Format: `"UNIQUE-ID | Title (Author1, Author2 et al., Year)"`
   - Includes title (truncated to 100 chars), first 3 authors, and year
   - Why: Provides human-readable citation text when PubMed/DOI unavailable - users can copy and search manually and later lookup might be possible


---

## License

- The license for the data is available at https://plantcyc.org/?webform=license-agreement or at [LICENSE.txt](LICENSE.txt). This includes the data from PlantCyc transformed in GPML format and available in this repository [/biocyc_pathways](./biocyc_pathways_20251217115329/)
- The license for the code is available at [GNU AGPL 3.0](https://www.gnu.org/licenses/agpl-3.0)


## References

- BioCyc Database Collection: [https://biocyc.org/BioCycUserGuide.shtml](https://biocyc.org/BioCycUserGuide.shtml)
- GPML2021 Specification: [GPML Documentation](https://pathvisio.org/documentation/GPML)
- WikiPathways (Agrawal et al. 2024): [https://doi.org/10.1093/nar/gkad960](https://doi.org/10.1093/nar/gkad960)
- PathVisio: [pathvisio.org/](https://pathvisio.org/)
- Plant Metabolic Network 16 (Hawkins et al. 2025) : [https://doi.org/10.1093/nar/gkae991](https://doi.org/10.1093/nar/gkae991)

## Contact / Maintainer 

Elena Del Pup, Wageningen University & Research 
elena.delpup@wur.nl 

## Use of AI Assistance
Development of this code involved  assistance from AI tools (Anthropic Claude Code, Opus, and Sonnet 4.5). These tools were used to speed up development, provide debugging suggestions, and help structure portions of the code. All AI-generated content was thoroughly reviewed, verified, and refined by the author.
