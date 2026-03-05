# Direct-Evolution-Monitoring-Web-Portal
MSc Bioinformatics Software Development Group Project 2026

## Variant Analysis and Activity Scoring
#### *Final script - starryeyed6.md/starryeyed6.py*

Python pipelines for mutation detection and activity scoring in directed evolution variant libraries.

Theese pipelines process plasmid sequence data and experimental measurements to identify mutations, classify their effects, and compute normalised activity scores relative to wild-type controls.

### Overview

Directed evolution experiments often generate lareg libarraies of plasmid variants and this script provides a copmputational workflow that:
- identifies coding sequences within circular plasmids
- translates variant genes into proteins
- aligns variant proteins against a wild-types reference
- detects substitutions, insertions, deletions, and truncations
- computes normalised activity scores based on DNA and protein quantification

The resulting dataset links sequence variation with experimental activity measurements, enabling systematic analysis of variant performance. 

### Pipeline Workflow

Input dataset -> Quality Control -> ORF Identification -> DNA → Protein Translation -> Protein Alignment (WT vs Variant) -> Mutation Annotation -> Variant Analysis -> Activity Scoring -> Annotated Output Dataset

### Key Features

- Handles circular plasmid sequences
- Robust ORF detection using approximate seed matching
- Protein-level global sequnce alignment
- Automatic mutation annotation
- Codon-level mutation classification
- Normalised activity scoring relative to the wild-type

### Input Data Format 

The pipeline accepts **TSV** or **JSON** datasets containing plasmid variants and experimental measurements

**Required Columns**

| Column                          | Description                                |
| ------------------------------- | ------------------------------------------ |
| `Plasmid_Variant_Index`         | Unique variant identifier                  |
| `Parent_Plasmid_Variant`        | Parent plasmid ID                          |
| `Directed_Evolution_Generation` | Evolution generation number                |
| `Assembled_DNA_Sequence`        | Circular plasmid DNA sequence              |
| `DNA_Quantification_fg`         | DNA output measurement                     |
| `Protein_Quantification_pg`     | Protein abundance measurement              |
| `Control`                       | Boolean flag indicating wild-type controls |


### Mutation Detection 

The pipeline identifies mutations by:
- Extracting the coding sequence (CDS) from each plasmid
- Translating CDS into protein sequence
- Performing global alignment with the wild-type protein
- Annotating differences using standard mutation notation

  *Example mutations annotations:*
- *T599S*
- *A600del*
- *ins601G*

Mutations are classifies as: 
- synonymous
- nonsynonymous
- insertions
- deletions
- truncating mutations

### Activity Scoring

Variant activity is quantifies using normalised DNA output per protein

**Raw activity**
Raw Activity = DNA_Quantification / Protein_Quantification

**Normalised activity score**
Activity Score =

    Raw Activity_variant
    --------------------
      Raw Activity_WT

*Interpretation*

| Score | Interpretation     |
| ----- | ------------------ |
| ≈ 1   | Wild-type activity |
| > 1   | Improved variant   |
| < 1   | Reduced activity   |

Two normalisations modes are supported:
- **global baseline** (all WT controls)
- **generation-specific** baseline

### Output 

The pipeline generates annotated datasets containing both experimental measurements and computed mutation data.

*Output files:*
- variant_analysis_output.csv
- variant_analysis_output.json
- variant_analysis.log
- activity_scores.csv
- activity_scores.json

Additional columns include:

| Column           | Description                        |
| ---------------- | ---------------------------------- |
| `protein_seq`    | Translated protein sequence        |
| `mutations`      | Mutation annotations               |
| `mutation_count` | Total mutation count               |
| `synonymous`     | Synonymous mutation count          |
| `nonsynonymous`  | Nonsynonymous mutation count       |
| `insertions`     | Number of insertions               |
| `deletions`      | Number of deletions                |
| `truncating`     | Premature stop codon indicator     |
| `raw_activity`   | DNA / protein activity metric      |
| `activity_score` | Normalized activity relative to WT |


### Installation dependencies using pip:

pip install -r requirements.txt

**requirements.txt*
- pandas==2.3.3
- biopython==1.86

### Software Environment

*Tested with:*
- Python 3.12.12
- pandas 2.3.3
- biopython 1.86

*Standard library modules used:*
- logging
- sys

Sequence alignment uses Biopython's PairwiseAligner API

### Example Usage

Run the analysis pipeline:
variant_analysis("input_dataset.tsv", best_orf)

Then compute activity scores:
activity_scoring(df, baseline_mode="global")


**input dataset can be either .tsv or .json*

**'df' in activity_scoring() should be output of variant_analysis() pipeline*

**'baseline_mode' can be either 'global' or 'generation'*
