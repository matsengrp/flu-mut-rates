# flu-syn-rates

## Overview

This pipeline calculates synonymous and non-synonymous mutation rates across influenza virus sequences by analyzing phylogenetic trees and sequence data. It processes multiple influenza subtypes (e.g., H1, H3, N1, N2), genome segments (e.g., HA, NA, PB2), and host groups (human, avian) as specified in the configuration.

## Setup

### Directory Structure

flu-syn-rates/
├── Snakefile             # Main workflow definition
├── config.yaml           # Configuration file
├── scripts/
│   ├── make_coding_sites.py        # Generate coding sites file
│   ├── make_count_dfs.py           # Count mutations along phylogenetic trees
│   ├── align_proteins.py           # Align proteins across subtypes
│   ├── rates_model.py              # Fit neutral models to mutation rates
│   ├── augment_expected_rates.py   # Add CG/GC to expected rates table
│   ├── MATWrapper.py               # Interface to MAT tree format
│   ├── ExpectedCalc.py             # Calculate expected mutation counts
│   └── basic_plots.py              # Plotting utilities
├── notebooks/
│   ├── compute_rates.ipynb              # Calculate mutation rates
│   └── analyze_genome_wide_rates.ipynb  # Visualize and analyze rates
├── logs/                 # Log files for pipeline runs
├── data/                 # Input data (organized by segment, subtype, and host)
│   ├── HA/
│   │   ├── H1/
│   │   │   ├── curated_root.fasta
│   │   │   ├── curated_reference.gff
│   │   │   ├── curated_reference.gtf
│   │   │   ├── final_tree.pb.gz          # Global tree (all hosts)
│   │   │   └── host_specific_trees/
│   │   │       ├── human_tree.pb.gz
│   │   │       ├── avian_tree.pb.gz
│   │   │       └── swine_tree.pb.gz
│   │   └── H3/
│   │       └── ...
│   └── NA/
│       └── ...
└── results/              # Output directory

### Configuration

Edit `config.yaml` to specify:
- HA and NA subtypes to analyze
- Genome segments to analyze
- Host groups to analyze
- Data directory locations

Example configuration:

```yaml
ha_subtypes:
  - "H1"
  - "H3"
  - "H5"

na_subtypes:
  - "N1"
  - "N2"

segments:
  - "HA"
  - "NA"
  - "PB2"
  # ... other segments

host_groups:
  - "human"
  - "avian"

data_dir: "../flu-usher/results"
output_dir: "results"
```

## Pipeline Steps

### Step 1: Coding Sites Identification

For each segment and subtype combination, the pipeline:
1. Reads a reference sequence (FASTA) and gene annotation (GFF)
2. Identifies coding and non-coding regions
3. Maps each nucleotide site to its codon position (1, 2, 3 or "noncoding")
4. Generates a CSV file with site-specific information including gene, codon site, and codon position

Note: Coding sites are shared across all host groups for a given segment/subtype combination.

### Step 2: Mutation Counting

The pipeline performs two types of mutation counting:

#### Global Tree Analysis (All Hosts Combined)
For each segment and subtype combination:
1. Takes the coding sites file from Step 1
2. Reads the global phylogenetic tree (`final_tree.pb.gz`), reference sequence, and GTF annotation
3. Traverses the tree to count mutations at each branch across all hosts
4. Classifies mutations as synonymous or non-synonymous
5. Outputs two files at the segment/subtype level:
   - `mutation_counts.csv` - Summary of mutations by site and type
   - `parent_child_pairs.csv` - Detailed mutation records for each branch

#### Host-Specific Tree Analysis
For each segment, subtype, and host combination:
1. Takes the coding sites file from Step 1
2. Reads the host-specific phylogenetic tree (`host_specific_trees/{host}_tree.pb.gz`), reference sequence, and GTF annotation
3. Traverses the tree to count mutations at each branch for that specific host
4. Classifies mutations as synonymous or non-synonymous
5. Outputs one file per host group:
   - `mutation_counts.csv` - Summary of mutations by site and type

### Step 3: Protein Alignment

For HA and NA segments only:
1. Aligns protein sequences across all subtypes for cross-subtype comparisons
2. Uses MUSCLE for multiple sequence alignment
3. Outputs aligned protein files to `results/aligned_proteins/{segment}/`

### Step 4: Compute Mutation Rates

Aggregates all mutation count data and calculates rates:
1. Combines mutation counts from all segments, subtypes, and hosts
2. Computes genome-wide mutation rates by mutation type and class (synonymous, nonsynonymous, nonsense)
3. Calculates segment-wide rates for each segment
4. Computes motif-level rates (3-mer sequence context)
5. Determines evolutionary opportunity thresholds
6. Calculates site-specific mutation rates

Outputs six CSV files at the results root:
- `counts.csv` - Aggregated mutation counts
- `genome_wide_rates.csv` - Genome-wide mutation rates
- `segment_wide_rates.csv` - Segment-specific mutation rates
- `motif_level_genome_wide_rates.csv` - Context-dependent rates
- `evo_opp_thresholds.csv` - Evolutionary opportunity filters
- `site_specific_mutation_rates.csv` - Per-site mutation rates

### Step 5: Analyze Genome-Wide Rates

Executes an analysis notebook that:
1. Visualizes synonymous mutation rates across mutation types
2. Compares rates between hosts (human vs avian)
3. Analyzes nonsynonymous and nonsense mutation rates relative to synonymous rates
4. Compares influenza rates with SARS-CoV-2 and HIV
5. Examines symmetry between complement mutation types

### Step 6: Fit Neutral Models

Fits log-linear models to predict synonymous mutation rates:
1. **Base model**: No context factors (mutation type only)
2. **Local context model**: Includes 3-mer sequence motif
3. **Full model**: Includes both local context (motif) and global context (segment)

Each model outputs:
- `expected_rates_by_predictor.csv` - Predicted rates for each combination of factors
- `model_performance.csv` - Model fit statistics (MSE, variance, R²)

Note: CG and GC mutation types are excluded from model fitting due to insufficient data.

### Step 7: Augment Expected Rates

Creates a complete expected rates table:
1. Loads fitted rates from the full model (local + global context)
2. Adds missing CG and GC mutation types using segment-specific empirical rates
3. Generates all 16 possible 3-mer motifs for CG and GC
4. Validates that all 12 mutation types have balanced entry counts
5. Outputs `expected_rates.csv` with complete coverage of all mutation types

## Running the Pipeline

Execute the full pipeline with:

```bash
snakemake --cores <N>
```

Or run specific targets:

```bash
# Process global tree for a specific segment/subtype
snakemake --cores 8 results/HA/H1/mutation_counts.csv results/HA/H1/parent_child_pairs.csv

# Process host-specific tree for a specific segment/subtype/host
snakemake --cores 8 results/HA/H1/human/mutation_counts.csv

# Process both global and host-specific for a segment/subtype
snakemake --cores 8 results/HA/H1/mutation_counts.csv results/HA/H1/human/mutation_counts.csv results/HA/H1/avian/mutation_counts.csv
```

## Output Files

The pipeline generates the following output files:

### Per-Segment/Subtype Outputs

1. **Coding Sites File**: `results/{segment}/{subtype}/coding_sites.csv`
   - Maps each nucleotide position to its coding properties
   - Shared across all analyses for a given segment/subtype combination

2. **Global Tree Outputs**: `results/{segment}/{subtype}/`
   - `mutation_counts.csv` - Aggregated mutation counts across all hosts
   - `parent_child_pairs.csv` - Detailed branch-level mutation information

3. **Host-Specific Tree Outputs**: `results/{segment}/{subtype}/{host}/`
   - `mutation_counts.csv` - Aggregated mutation counts for the specific host group

4. **Aligned Proteins** (HA and NA only): `results/aligned_proteins/{segment}/`
   - Cross-subtype protein alignments

### Genome-Wide Analysis Outputs

Located in the `results/` root directory:

5. **Mutation Rate Files**:
   - `counts.csv` - Combined mutation counts from all segments and subtypes
   - `genome_wide_rates.csv` - Mutation rates by type and class (synonymous, nonsynonymous, nonsense)
   - `segment_wide_rates.csv` - Segment-specific mutation rates
   - `motif_level_genome_wide_rates.csv` - Context-dependent rates (3-mer motifs)
   - `evo_opp_thresholds.csv` - Evolutionary opportunity thresholds for filtering
   - `site_specific_mutation_rates.csv` - Per-site mutation rates

6. **Neutral Model Outputs**: `results/neutral_model/`
   - `base/` - Model with mutation type only
     - `expected_rates_by_predictor.csv`
     - `model_performance.csv`
   - `local_context/` - Model with mutation type + 3-mer motif
     - `expected_rates_by_predictor.csv`
     - `model_performance.csv`
   - `local_context+global_context/` - Model with mutation type + motif + segment
     - `expected_rates_by_predictor.csv`
     - `model_performance.csv`

7. **Complete Expected Rates**: `results/expected_rates.csv`
   - Full model predictions augmented with CG and GC empirical rates
   - Contains all 12 mutation types × 16 motifs × 8 segments (1536 rows)

### Output Structure Example

```
results/
├── HA/
│   ├── H1/
│   │   ├── coding_sites.csv
│   │   ├── mutation_counts.csv
│   │   ├── parent_child_pairs.csv
│   │   ├── human/
│   │   │   └── mutation_counts.csv
│   │   ├── avian/
│   │   │   └── mutation_counts.csv
│   │   └── swine/
│   │       └── mutation_counts.csv
│   └── H3/
│       └── ...
├── aligned_proteins/
│   ├── HA/
│   └── NA/
├── counts.csv
├── genome_wide_rates.csv
├── segment_wide_rates.csv
├── motif_level_genome_wide_rates.csv
├── evo_opp_thresholds.csv
├── site_specific_mutation_rates.csv
├── expected_rates.csv
└── neutral_model/
    ├── base/
    │   ├── expected_rates_by_predictor.csv
    │   └── model_performance.csv
    ├── local_context/
    │   ├── expected_rates_by_predictor.csv
    │   └── model_performance.csv
    └── local_context+global_context/
        ├── expected_rates_by_predictor.csv
        └── model_performance.csv
```

## Requirements

### Software Dependencies

- Python 3.10+
- Snakemake
- Jupyter (for notebook execution)
- MUSCLE (for protein alignment)
- Standard scientific Python packages (pandas, numpy, scipy, matplotlib, seaborn, biopython)

### Environment Setup

Create the conda environment from the provided file:

```bash
conda env create -f environment.yml
conda activate flu-syn-rates
```

The pipeline will automatically execute Jupyter notebooks as part of the workflow using `jupyter nbconvert`.

