# flu-syn-rates

## Overview

This pipeline calculates synonymous and non-synonymous mutation rates across influenza virus sequences by analyzing phylogenetic trees and sequence data. It processes multiple influenza subtypes (e.g., H1, H3, N1, N2), genome segments (e.g., HA, NA, PB2), and host groups (human, avian) as specified in the configuration.

## Setup

### Directory Structure

flu-syn-rates/
├── Snakefile             # Main workflow definition
├── config/
│   └── config.yaml       # Configuration file
├── scripts/
│   ├── make_coding_sites.py   # Script to generate coding sites file
│   └── make_count_dfs.py      # Script to count mutations along tree
├── data/                 # Input data (organized by segment, subtype, and host)
│   ├── HA/
│   │   ├── H1/
│   │   │   ├── curated_root.fasta
│   │   │   ├── curated_reference.gff
│   │   │   ├── curated_reference.gtf
│   │   │   └── host_specific_trees/
│   │   │       ├── human_tree.pb.gz
│   │   │       └── avian_tree.pb.gz
│   │   └── H3/
│   │       └── ...
│   └── NA/
│       └── ...
└── results/              # Output directory

### Configuration

Edit `config/config.yaml` to specify:
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

For each segment, subtype, and host combination, the pipeline:
1. Takes the coding sites file from Step 1
2. Reads the host-specific phylogenetic tree, reference sequence, and GTF annotation
3. Traverses the tree to count mutations at each branch
4. Classifies mutations as synonymous or non-synonymous
5. Outputs two files per host group:
   - `mutation_counts.csv` - Summary of mutations by site and type
   - `parent_child_pairs.csv` - Detailed mutation records for each branch

## Running the Pipeline

Execute the full pipeline with:

```bash
snakemake --cores <N>
```

Or run specific targets:

```bash
# Process a specific segment/subtype/host combination
snakemake --cores 8 results/HA/H1/human/mutation_counts.csv

# Process all host groups for a specific segment/subtype
snakemake --cores 8 results/HA/H1/human/mutation_counts.csv results/HA/H1/avian/mutation_counts.csv
```

## Output Files

The pipeline generates the following output files:

1. **Coding Sites File**: `results/{segment}/{subtype}/coding_sites.csv`
   - Maps each nucleotide position to its coding properties
   - Shared across all host groups for a given segment/subtype combination

2. **Host-Specific Mutation Count Files**: `results/{segment}/{subtype}/{host}/`
   - `mutation_counts.csv` - Aggregated mutation counts for the host group
   - `parent_child_pairs.csv` - Detailed mutation information for the host group

3. **Aligned Proteins** (HA and NA only): `results/aligned_proteins/{segment}/`
   - Cross-subtype protein alignments

Example output structure:
```
results/
├── HA/
│   ├── H1/
│   │   ├── coding_sites.csv          # Shared
│   │   ├── human/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   └── avian/
│   │       ├── mutation_counts.csv
│   │       └── parent_child_pairs.csv
│   └── H3/
│       └── ...
└── aligned_proteins/
    └── HA/
```

## Requirements

