# flu-syn-rates

## TODO
* Analyze fitness effects
   * compare to DMS data
   * visualize fitness effects

## Overview

This pipeline calculates synonymous and non-synonymous mutation rates across influenza virus sequences by analyzing phylogenetic trees and sequence data. It processes multiple influenza subtypes (e.g., H1N1, H3N2) and genome segments (e.g., HA, NA) as specified in the configuration.

## Setup

### Directory Structure

flu-syn-rates/
├── Snakefile             # Main workflow definition
├── config/
│   └── config.yaml       # Configuration file
├── scripts/
│   ├── make_coding_sites.py   # Script to generate coding sites file
│   └── make_count_dfs.py      # Script to count mutations along tree
├── data/                 # Input data (organized by subtype and segment)
│   ├── H1N1/
│   │   ├── HA/
│   │   │   ├── reference.fasta
│   │   │   ├── annotation.gff
│   │   │   ├── annotation.gtf
│   │   │   └── tree.nwk
│   │   └── NA/
│   │       └── ...
│   └── H3N2/
│       └── ...
└── results/              # Output directory

### Configuration

Edit `config/config.yaml` to specify:
- Subtypes and segments to analyze
- Data directory locations
- File naming patterns

Example configuration:

subtypes:
  - "H1N1"
  - "H3N2"

segments:
  - "HA"
  - "NA"

data_dir: "data/"
output_dir: "results/"

reference_format: "{data_dir}/{subtype}/{segment}/reference.fasta"
gff_format: "{data_dir}/{subtype}/{segment}/annotation.gff"
gtf_format: "{data_dir}/{subtype}/{segment}/annotation.gtf"
tree_format: "{data_dir}/{subtype}/{segment}/tree.nwk"

## Pipeline Steps

### Step 1: Coding Sites Identification

For each subtype and segment combination, the pipeline:
1. Reads a reference sequence (FASTA) and gene annotation (GFF)
2. Identifies coding and non-coding regions
3. Maps each nucleotide site to its codon position (1, 2, 3 or "noncoding")
4. Generates a CSV file with site-specific information including gene, codon site, and codon position

### Step 2: Mutation Counting

For each subtype and segment combination, the pipeline:
1. Takes the coding sites file from Step 1
2. Reads the phylogenetic tree, reference sequence, and GTF annotation
3. Traverses the tree to count mutations at each branch
4. Classifies mutations as synonymous or non-synonymous
5. Outputs two files:
   - `mutation_counts.csv` - Summary of mutations by site and type
   - `parent_child_pairs.csv` - Detailed mutation records for each branch

## Running the Pipeline

Execute the pipeline with:

snakemake --use-conda --cores <N>

## Output Files

The pipeline generates two main output files for each subtype/segment:

1. **Coding Sites File**: `results/{subtype}/{segment}/coding_sites.csv`
   - Maps each nucleotide position to its coding properties

2. **Mutation Count Files**:
   - `results/{subtype}/{segment}/mutation_counts.csv` - Aggregated mutation counts
   - `results/{subtype}/{segment}/parent_child_pairs.csv` - Detailed mutation information

## Requirements

