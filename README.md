# flu-syn-rates

## Overview

This pipeline calculates synonymous and non-synonymous mutation rates across influenza virus sequences by analyzing phylogenetic trees and sequence data. It processes multiple influenza subtypes (e.g., H1, H3, N1, N2), genome segments (e.g., HA, NA, PB2), and host groups (human, avian) as specified in the configuration.

## Setup

### Directory Structure

```
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
```

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
2. Adds missing CG and GC mutation types using empirical rates: motif-specific genome-wide rates are used where available (12 of 16 motifs per type), with segment-wide rates as a fallback for the remaining motifs
3. Generates all 16 possible 3-mer motifs for CG and GC
4. Validates that all 12 mutation types have balanced entry counts
5. Outputs `expected_rates.csv` with complete coverage of all mutation types

### Step 8: Compute Fitness Effects

Computes fitness effects of synonymous and amino acid mutations by comparing observed mutation counts to neutral model expectations:
1. Joins `counts.csv` with `expected_rates.csv` to compute expected counts per site under the neutral model
2. Aggregates synonymous mutations by site and wildtype nucleotide to compute per-site synonymous fitness effects
3. Aggregates all nucleotide mutations by resulting amino acid change to compute per-amino-acid-mutation fitness effects
4. Fitness effects are estimated as log((actual_count + 0.5) / (expected_count + 0.5))

Outputs three CSV files at the results root:
- `actual_expected.csv` - Per-mutation counts joined with neutral model expected counts
- `sitewise_synonymous_fitness_effects.csv` - Per-site synonymous fitness effects
- `aa_fitness_effects.csv` - Per-amino-acid-mutation fitness effects

### Step 9: Process DMS Data

Processes raw deep mutational scanning (DMS) data from external sources into standardized formats for comparison with fitness effects. Each experiment is handled by a separate notebook:

**`process_dms_data_yu_ha.ipynb`** (Yu et al., HA):
1. Aligns the DMS experiment HA sequence to the H3 tree reference sequence using MUSCLE to establish site numbering correspondence
2. Merges DMS phenotype measurements with the numbering map
3. Outputs `results/dms_data/Yu_HA/processed_dms_data.csv` — processed HA DMS data with tree reference site numbering

**`process_dms_data_bloom_np.ipynb`** (Bloom et al., NP):
1. Verifies the DMS sequence matches the NP tree reference (Aichi 1968)
2. Computes log-ratio mutation effects (log(preference / wt_preference))
3. Outputs `results/dms_data/Bloom_NP/processed_dms_data.csv` — NP DMS data with log-ratio fitness effects

**`process_dms_data_soh_pb2.ipynb`** (Soh et al., PB2):
1. Aligns the DMS sequence to the PB2 reference using MUSCLE
2. Verifies the alignment has no gaps and reports percent identity (QC only; raw data is used directly by Step 10)

**`process_dms_data_wang_na.ipynb`** (Wang et al., NA):
1. Reads NA DMS data and compares the DMS sequence to the N1 tree reference (sequence comparison / exploration)

### Step 10: Analyze Fitness Effects

Executes an analysis notebook that:
1. Plots distributions of fitness effects by mutation class across all genes
2. Examines synonymous fitness effects across genome sites
3. Compares evolutionary fitness effects to experimentally measured DMS effects for HA (Yu et al.), NP (Bloom et al.), and PB2 (Soh et al.)

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
   - Columns:
     - `site` — nucleotide position in the reference sequence (0-indexed)
     - `codon_position` — position within the codon (1, 2, 3, or `"noncoding"`)
     - `codon_site` — index of the codon within the gene
     - `gene` — gene name the site belongs to (or `"noncoding"`)

2. **Global Tree Outputs**: `results/{segment}/{subtype}/`
   - `mutation_counts.csv` - Aggregated mutation counts across all hosts
     - Columns:
       - `site` — nucleotide position in the reference sequence (1-indexed)
       - `nt_mut` — nucleotide mutation string (e.g. `"A1C"`: wildtype + position + mutant)
       - `wt_nt` — wildtype nucleotide
       - `mut_nt` — mutant nucleotide
       - `gene` — gene containing the site
       - `codon_position` — position within the codon (1, 2, or 3)
       - `codon_site` — index of the codon within the gene (1-indexed)
       - `wt_codon` — wildtype codon sequence
       - `mut_codon` — mutant codon sequence
       - `wt_aa` — wildtype amino acid (single-letter code)
       - `mut_aa` — mutant amino acid (single-letter code)
       - `aa_mut` — amino acid mutation string (e.g. `"M1L"`: wildtype AA + codon site + mutant AA)
       - `parent_motif` — 3-mer sequence context of the parent node at the mutated site
       - `actual_count` — number of times this mutation was observed across tree branches
       - `branch_length` — total branch length over which the mutation opportunity was counted
       - `mut_class` — mutation class (`"synonymous"`, `"nonsynonymous"`, `"nonsense"`, or `"noncoding"`)
       - `mut_type` — two-character mutation type (e.g. `"AC"` for A→C)
   - `parent_child_pairs.csv` - Detailed branch-level mutation information
     - Columns:
       - `parent_name` — identifier of the parent node in the phylogenetic tree
       - `child_name` — identifier of the child node in the phylogenetic tree
       - `parent` — full nucleotide sequence reconstructed at the parent node
       - `child` — full nucleotide sequence reconstructed at the child node
       - `branch_length` — branch length between parent and child nodes

3. **Host-Specific Tree Outputs**: `results/{segment}/{subtype}/{host}/`
   - `mutation_counts.csv` - Aggregated mutation counts for the specific host group (same columns as global `mutation_counts.csv`)

4. **Aligned Proteins** (HA and NA only): `results/aligned_proteins/{segment}/`
   - Cross-subtype protein alignments

### Genome-Wide Analysis Outputs

Located in the `results/` root directory:

5. **Mutation Rate Files**:
   - `counts.csv` - Combined mutation counts from all segments and subtypes
     - Columns: same as `mutation_counts.csv` except `parent_motif` is renamed to `motif`, plus:
       - `motif` — 3-mer sequence context of the mutation (centered on the mutated site)
       - `subtype` — influenza subtype (e.g. `"H1"`, `"N2"`, or `"all"`)
       - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
       - `segment_subtype` — combined segment and subtype label (e.g. `"HA_H1"`)
       - `segment_length` — length of the genome segment in nucleotides
       - `host` — host group (`"human"`, `"avian"`, `"swine"`, or `"all"`)
       - `evo_opp` — evolutionary opportunity
       - `rate` — per-site mutation rate
   - `genome_wide_rates.csv` - Mutation rates by type and class (synonymous, nonsynonymous, nonsense)
     - Columns:
       - `mut_type` — two-character mutation type (e.g. `"AC"` for A→C)
       - `mut_class` — mutation class (`"synonymous"`, `"nonsynonymous"`, or `"nonsense"`)
       - `host` — host group (`"human"`, `"avian"`, `"swine"`, or `"all"`)
       - `actual_count` — total number of observed mutations of this type/class/host
       - `evo_opp` — evolutionary opportunity
       - `rate` — mutation rate
       - `syn_rate` — synonymous mutation rate for the same host and mutation type
       - `rel_rate` — rate relative to the synonymous rate (`rate / syn_rate`)
   - `segment_wide_rates.csv` - Segment-specific mutation rates
     - Columns: same as `genome_wide_rates.csv` minus `syn_rate` and `rel_rate`, plus:
       - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
   - `motif_level_genome_wide_rates.csv` - Context-dependent rates (3-mer motifs)
     - Columns:
       - `mut_type` — two-character mutation type (e.g. `"AC"` for A→C)
       - `motif` — 3-mer sequence context (centered on the mutated site)
       - `host` — host group (`"human"`, `"avian"`, `"swine"`, or `"all"`)
       - `actual_count` — total number of observed mutations for this type/motif/host
       - `evo_opp` — evolutionary opportunity
       - `rate` — mutation rate
   - `evo_opp_thresholds.csv` - Evolutionary opportunity thresholds for filtering
     - Columns: same as `genome_wide_rates.csv` minus `rel_rate`, plus:
       - `evo_opp_threshold` — minimum evolutionary opportunity a site must have to be included in rate calculations for this mutation type/class/host
   - `site_specific_mutation_rates.csv` - Per-site mutation rates
     - Columns: same as `counts.csv`, plus:
       - `evo_opp_threshold` — minimum evolutionary opportunity threshold for this mutation type/class/host (from `evo_opp_thresholds.csv`)
       - `unclipped_rate` — raw per-site mutation rate before applying the floor
       - `min_rate` — floor rate applied to sites below the evolutionary opportunity threshold (set to the synonymous genome-wide rate)

6. **Neutral Model Outputs**: `results/neutral_model/`
   - `base/` - Model with mutation type only
     - `expected_rates_by_predictor.csv`
       - Columns:
         - `mut_type` — two-character mutation type (e.g. `"AC"` for A→C)
         - `predicted_rate` — model-predicted synonymous mutation rate
     - `model_performance.csv`
       - Columns (same for all three models):
         - `mut_type` — two-character mutation type (e.g. `"AC"` for A→C)
         - `mse` — mean squared error of the model fit
         - `var` — variance of the observed rates (baseline MSE with no predictors)
         - `r2` — R² of the model fit
   - `local_context/` - Model with mutation type + 3-mer motif
     - `expected_rates_by_predictor.csv`
       - Columns: same as `base/expected_rates_by_predictor.csv`, plus:
         - `motif` — 3-mer sequence context
     - `model_performance.csv`
   - `local_context+global_context/` - Model with mutation type + motif + segment
     - `expected_rates_by_predictor.csv`
       - Columns: same as `local_context/expected_rates_by_predictor.csv`, plus:
         - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
     - `model_performance.csv`

7. **Complete Expected Rates**: `results/expected_rates.csv`
   - Full model predictions augmented with CG and GC empirical rates
   - Contains all 12 mutation types × 16 motifs × 8 segments (1536 rows)
   - Columns: same as `local_context+global_context/expected_rates_by_predictor.csv` (`mut_type`, `segment`, `motif`, `predicted_rate`), but with complete coverage of all 12 mutation types; CG and GC rates are filled in from the augmentation step (Step 7) rather than model predictions

8. **Fitness Effect Files**:
   - `actual_expected.csv` - Per-mutation counts joined with neutral model expected counts
     - Columns: same as `counts.csv`, plus:
       - `predicted_rate` — neutral model predicted rate for this mutation type/motif/segment (from `expected_rates.csv`)
       - `expected_count` — expected number of observations under the neutral model (`predicted_rate × evo_opp`)
   - `sitewise_synonymous_fitness_effects.csv` - Per-site synonymous fitness effects
     - Columns:
       - `host` — host group (`"human"`, `"avian"`, `"swine"`, or `"all"`)
       - `subtype` — influenza subtype (e.g. `"H1"`, `"N2"`, or `"all"`)
       - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
       - `gene` — gene containing the site
       - `site` — nucleotide position in the reference sequence (1-indexed)
       - `codon_site` — index of the codon within the gene (1-indexed)
       - `wt_nt` — wildtype nucleotide at this site
       - `actual_count` — total observed synonymous mutations away from `wt_nt` at this site
       - `expected_count` — total expected synonymous mutations under the neutral model
       - `delta_fitness` — estimated fitness effect: log((actual_count + 0.5) / (expected_count + 0.5))
   - `aa_fitness_effects.csv` - Per-amino-acid-mutation fitness effects
     - Columns:
       - `host` — host group (`"human"`, `"avian"`, `"swine"`, or `"all"`)
       - `subtype` — influenza subtype (e.g. `"H1"`, `"N2"`, or `"all"`)
       - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
       - `gene` — gene containing the site
       - `codon_site` — index of the codon within the gene (1-indexed)
       - `wt_aa` — wildtype amino acid (single-letter code)
       - `mut_aa` — mutant amino acid (single-letter code)
       - `aa_mut` — amino acid mutation string (e.g. `"M1L"`)
       - `mut_class` — mutation class (`"synonymous"`, `"nonsynonymous"`, or `"nonsense"`)
       - `actual_count` — total observed nucleotide mutations resulting in this amino acid change
       - `expected_count` — total expected mutations under the neutral model
       - `delta_fitness` — estimated fitness effect: log((actual_count + 0.5) / (expected_count + 0.5))

9. **Processed DMS Data**: `results/dms_data/`
   - `Yu_HA/processed_dms_data.csv` - HA DMS phenotypes (Yu et al.) with tree reference site numbering
   - `Bloom_NP/processed_dms_data.csv` - NP DMS preferences (Bloom et al.) with log-ratio fitness effects
   - `.process_dms_data_soh_pb2.done` - Sentinel file marking completion of Soh et al. PB2 alignment QC
   - `.process_dms_data_wang_na.done` - Sentinel file marking completion of Wang et al. NA sequence comparison

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
├── actual_expected.csv
├── sitewise_synonymous_fitness_effects.csv
├── aa_fitness_effects.csv
├── dms_data/
│   ├── Yu_HA/
│   │   └── processed_dms_data.csv
│   └── Bloom_NP/
│       └── processed_dms_data.csv
├── .process_dms_data_soh_pb2.done
├── .process_dms_data_wang_na.done
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

