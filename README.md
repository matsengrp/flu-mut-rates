# flu-mut-rates

## Overview

This pipeline calculates synonymous and non-synonymous mutation rates across influenza virus sequences by analyzing phylogenetic trees and sequence data. It processes multiple influenza subtypes (e.g., H1, H3, N1, N2), genome segments (e.g., HA, NA, PB2), and tree subsets (host groups, geographic regions, temporal periods) as specified in the configuration.

## Setup

### Directory Structure

```
flu-mut-rates/
├── Snakefile             # Main workflow definition
├── config.yaml           # Configuration file
├── scripts/
│   ├── make_coding_sites.py        # Generate coding sites file
│   ├── make_count_dfs.py           # Count mutations along phylogenetic trees
│   ├── align_proteins.py           # Align proteins across subtypes
│   ├── rates_model.py              # Fit neutral models to mutation rates
│   ├── augment_expected_rates.py   # Add CG/GC to expected rates table
│   └── ExpectedCalc.py             # Calculate expected mutation counts
├── notebooks/
│   ├── compute_rates.ipynb                    # Calculate mutation rates
│   ├── analyze_genome_wide_rates.ipynb        # Visualize and analyze rates
│   ├── analyze_site_specific_rates.ipynb      # Site-specific rate analysis with SHAPE-MaP
│   ├── compute_fitness_effects.ipynb          # Compute per-mutation fitness effects
│   ├── process_shapemap_data.ipynb            # Process SHAPE-MaP reactivity data
│   ├── process_dms_data_yu_ha.ipynb           # Process Yu et al. HA DMS data
│   ├── process_dms_data_bloom_np.ipynb        # Process Bloom et al. NP DMS data
│   ├── process_dms_data_soh_pb2.ipynb         # Process Soh et al. PB2 DMS data
│   ├── process_dms_data_wang_na.ipynb         # Process Wang et al. NA DMS data
│   ├── process_dms_data_li_pb1.ipynb          # Process Li et al. PB1 DMS data
│   ├── process_dms_data_hom_m1.ipynb          # Process Hom et al. M1 DMS data
│   ├── process_dms_data_teo_nep.ipynb         # Process Teo et al. NEP DMS data
│   ├── analyze_fitness_effects.ipynb          # Compare fitness effects to DMS data
│   ├── summarize_filter_logs.ipynb            # Summarize mutation filter statistics
│   ├── compute_subset_rates.ipynb             # Compute rates for subset trees (host, geographic, temporal)
│   ├── compute_subset_fitness_effects.ipynb   # Compute AA fitness effects per subset
│   ├── analyze_subset_fitness_effects.ipynb   # Compare fitness effects between subsets
│   ├── check_subset_pcp_overlap.ipynb         # Check PCP overlap between subsets
│   └── compose_figures.ipynb                  # Assemble multi-panel figures from individual PNGs
├── logs/                 # Log files for pipeline runs
├── data/                 # Input data (organized by segment, subtype, and host)
│   ├── packaging_signal_boundaries.csv  # Packaging signal boundaries per segment (from Li et al. 2021)
│   ├── splice_site_boundaries.csv       # Canonical M2 and NEP splice-site coordinates (from Lamb & Lai 1980)
│   ├── HA/
│   │   ├── H1/
│   │   │   ├── curated_root.fasta
│   │   │   ├── curated_reference.gff
│   │   │   ├── curated_reference.gtf
│   │   │   ├── final_tree.pb.gz          # Global tree (all hosts)
│   │   │   ├── host_specific_trees/
│   │   │   │   ├── human_tree.pb.gz
│   │   │   │   └── avian_tree.pb.gz
│   │   │   ├── geographic_trees/
│   │   │   │   ├── north_america_tree.pb.gz
│   │   │   │   ├── europe_tree.pb.gz
│   │   │   │   └── asia_tree.pb.gz
│   │   │   └── temporal_trees/
│   │   │       ├── early_tree.pb.gz
│   │   │       └── late_tree.pb.gz
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
- Geographic groups to analyze
- Temporal groups to analyze
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

geographic_groups:
  - "north_america"
  - "europe"
  - "asia"

temporal_groups:
  - "early"
  - "late"

data_dir: "../flu-usher/results"
output_dir: "results"
```

### Input data files

Two CSV tables in `data/` define regions of the genome that are under additional non-neutral constraint on otherwise-synonymous mutations. They are consumed by `notebooks/compute_rates.ipynb` (to build the `exclude_from_mut_rate_analysis` flag used during rate aggregation and neutral-model fitting; see Step 4) and by `notebooks/analyze_fitness_effects.ipynb` (to overlay the affected regions on per-site fitness-effect plots; see Step 12).

**`data/packaging_signal_boundaries.csv`** — length of each segment's packaging signal at either vRNA terminus.
- **Columns:** `segment`, `end` (`3prime_vRNA` or `5prime_vRNA`), `nt` (packaging-signal length in nucleotides), `reference` (literature source).
- **Interpretation:** site coordinates in the pipeline are CDS-based and 1-indexed, running 5′→3′ along the positive-sense mRNA. The `3prime_vRNA` end corresponds to the CDS start, so a value `n3` flags CDS sites `[1, n3]`. The `5prime_vRNA` end corresponds to the CDS end, so a value `n5` flags CDS sites `[cds_len − n5 + 1, cds_len]`, where `cds_len` is the segment's CDS length (available from `results/{segment}/{subtype}/coding_sites.csv`). An `nt` value of 0 means no packaging signal is annotated at that terminus (e.g. PB2 at the 3′ vRNA end).

**`data/splice_site_boundaries.csv`** — canonical M2 (MP segment) and NEP (NS segment) splice junctions.
- **Columns:** `segment`, `product` (`M2` or `NEP`), `feature` (`5prime_splice_site` or `3prime_splice_site`), `position` (half-integer CDS coordinate), `reference`.
- **Interpretation:** positions are half-integers (e.g. `26.5`) because splicing cuts *between* two adjacent nucleotides rather than at a single site; the half-integer places the junction between sites `floor(position)` and `ceil(position)`. Rate aggregation (Step 4) flags a 10-nt window `|site − position| ≤ 5` around each junction (e.g. position `26.5` → sites 22–31); the fitness-effect overlays (Step 12) draw a vertical marker at the half-integer coordinate so the marker sits between the two flanking sites.

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
5. Outputs two files per host group:
   - `mutation_counts.csv` - Summary of mutations by site and type
   - `parent_child_pairs.csv` - Detailed mutation records for each branch

#### Geographic Tree Analysis
For each segment, subtype, and geographic group combination:
1. Takes the coding sites file from Step 1
2. Reads the geographic phylogenetic tree (`geographic_trees/{geo}_tree.pb.gz`), reference sequence, and GTF annotation
3. Traverses the tree to count mutations at each branch for that geographic region
4. Classifies mutations as synonymous or non-synonymous
5. Outputs two files per geographic group:
   - `mutation_counts.csv` - Summary of mutations by site and type
   - `parent_child_pairs.csv` - Detailed mutation records for each branch

#### Temporal Tree Analysis
For each segment, subtype, and temporal group combination:
1. Takes the coding sites file from Step 1
2. Reads the temporal phylogenetic tree (`temporal_trees/{temporal}_tree.pb.gz`), reference sequence, and GTF annotation
3. Traverses the tree to count mutations at each branch for that time period
4. Classifies mutations as synonymous or non-synonymous
5. Outputs two files per temporal group:
   - `mutation_counts.csv` - Summary of mutations by site and type
   - `parent_child_pairs.csv` - Detailed mutation records for each branch

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

Computes fitness effects of nucleotide and amino acid mutations by comparing observed mutation counts to neutral model expectations:
1. Joins `counts.csv` with `expected_rates.csv` to compute expected counts per site under the neutral model
2. Aggregates counts by nucleotide mutation to compute per-nucleotide-mutation fitness effects as log((actual_count + 0.5) / (expected_count + 0.5))
3. Computes per-site synonymous fitness effects as the mean of per-nucleotide-mutation `delta_fitness` values across synonymous mutations at each site
4. Aggregates all nucleotide mutations by resulting amino acid change to compute per-amino-acid-mutation fitness effects as log((actual_count + 0.5) / (expected_count + 0.5))

Outputs four CSV files at the results root:
- `actual_expected.csv` - Per-mutation counts joined with neutral model expected counts
- `nt_fitness_effects.csv` - Per-nucleotide-mutation fitness effects
- `sitewise_synonymous_fitness_effects.csv` - Per-site synonymous fitness effects (mean of per-nucleotide-mutation effects)
- `aa_fitness_effects.csv` - Per-amino-acid-mutation fitness effects

### Step 9: Process SHAPE-MaP Data

Processes SHAPE-MaP RNA secondary structure reactivity data from Dadonaite et al. 2019:
1. Reads per-nucleotide SHAPE reactivity values from the source Excel file
2. Aligns each segment's reactivity profile to the pipeline reference sequences using MUSCLE
3. Outputs `results/shapemap/all_data.csv` — SHAPE reactivities mapped to reference site coordinates for all segments (HA and NA excluded)

### Step 10: Process DMS Data

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
2. Verifies the alignment has no gaps and reports percent identity (QC only; raw data is used directly by Step 11)

**`process_dms_data_wang_na.ipynb`** (Wang et al., NA):
1. Aligns the DMS NA sequence to the N1 tree reference using MUSCLE to establish site numbering correspondence
2. Merges fitness measurements with the numbering map
3. Outputs `results/dms_data/Wang_NA/processed_dms_data.csv` — NA DMS data with tree reference site numbering

**`process_dms_data_li_pb1.ipynb`** (Li et al., PB1):
1. Aligns the DMS PB1 sequence to the tree reference using MUSCLE to establish site numbering correspondence
2. Merges fitness measurements with the numbering map
3. Outputs `results/dms_data/Li_PB1/processed_dms_data.csv` — PB1 DMS data with tree reference site numbering

**`process_dms_data_hom_m1.ipynb`** (Hom et al., M1):
1. Aligns the DMS M1 sequence to the tree reference using MUSCLE to establish site numbering correspondence
2. Converts amino acid preferences to log-ratio fitness effects
3. Outputs `results/dms_data/Hom_M1/processed_dms_data.csv` — M1 DMS data with tree reference site numbering

**`process_dms_data_teo_nep.ipynb`** (Teo et al., NEP):
1. Reconstructs the NEP reference protein from the NS segment reference by concatenating exon 1 (nt 1–30) and exon 2 (nt 503–838)
2. Aligns the DMS NEP sequence (PR8 strain) to the reference using MUSCLE; uses actual DMS site numbers for coordinate mapping to handle non-consecutive site numbering
3. Outputs `results/dms_data/Teo_NEP/processed_dms_data.csv` — NEP DMS data with tree reference site numbering

**`process_dms_data_chen_pa.ipynb`** (Chen et al., PA):
1. Parses the AA mutation column to extract wildtype AA, DMS site, and mutant AA; excludes indel and stop codon mutations
2. Averages the two no-drug fitness replicates (P1NO-1-fit, P1NO-2-fit) and groups by AA mutation (multiple nucleotide changes can produce the same AA change)
3. Computes log-scale DMS effects: `log(mean fitness without drug)`
4. Aligns the DMS sequence (first 240 AA of PA) to the PA tree reference using MUSCLE to establish site numbering correspondence
5. Outputs `results/dms_data/Chen_PA/processed_dms_data.csv` — PA DMS data with tree reference site numbering

### Step 11: Analyze Site-Specific Rates

Executes an analysis notebook that:
1. Visualizes per-site synonymous mutation rates across the genome
2. Examines the relationship between local sequence context (motif) and site-specific rates
3. Compares model performance (base, local context, full context) across segments
4. Integrates SHAPE-MaP RNA structure data to examine the relationship between secondary structure and mutation rates

### Step 12: Analyze Fitness Effects

Executes an analysis notebook that:
1. Plots distributions of fitness effects by mutation class across all genes
2. Examines per-site synonymous fitness effects across the genome, with overlays highlighting regions under non-coding selection: packaging signals at vRNA termini (from `data/packaging_signal_boundaries.csv`), the PA-X alternative ORF on segment PA, and canonical splice sites for M2 (MP) and NEP (NS) (from `data/splice_site_boundaries.csv`). Splice-site positions are stored as half-integer midpoints (e.g. 26.5) because splice events cut between two adjacent nucleotides rather than at a single site; the half-integer places the rendered vertical marker between the flanking sites on the figure
3. Compares evolutionary fitness effects to experimentally measured DMS effects for seven proteins: HA/H3 (Yu et al. 2025), NP (Bloom et al. 2014), M1 (Hom et al. 2019), NEP (Teo et al. 2024), PB1 (Li et al. 2023), PB2 (Soh et al. 2019), and NA/N1 (Wang et al. 2023)

### Step 13: Summarize Mutation Filter Logs

Executes a diagnostic notebook that parses the log files produced by the mutation-counting rules and consolidates filter statistics across all trees. For each tree (global and host-specific), the notebook reports:
- Total nodes and internal nodes
- Branch counts: total examined, passing, and filtered
- Mutation counts: total, passing, and filtered (broken down by reason: too many mutations, zero mutations, duplicate codon targets)

The notebook produces three figures:
1. Branch filter summary per tree: passing branch counts (log scale) and % passing, grouped by host
2. Mutation filter summary per tree: passing mutation counts (log scale) and % passing, grouped by host
3. Stacked bar breakdown of % mutations by filter reason, faceted by host

### Step 14: Compute Subset Rates

Aggregates mutation counts from all subset trees (host-specific, geographic, temporal):
1. Reads mutation counts from each subset directory
2. Labels each row with `subset` (group name) and `subset_type` (host/geographic/temporal)
3. Computes evolutionary opportunity and mutation rates
4. Outputs `results/subset_counts.csv`

### Step 15: Compute Subset Fitness Effects

Computes amino-acid-level fitness effects for each subset:
1. Reads `subset_counts.csv` and `expected_rates.csv` (from the global neutral model)
2. Merges to compute expected counts per mutation under the neutral model
3. Filters to mutations with at least 10 actual or expected counts
4. Groups by amino acid mutation and computes fitness effects as log((actual + 0.5) / (expected + 0.5))
5. Outputs `results/subset_aa_fitness_effects.csv`

### Step 16: Analyze Subset Fitness Effects

Compares fitness effects between subsets:
1. For each pair of subsets within the same type (e.g., human vs avian, early vs late, north_america vs europe):
   - Matches nonsynonymous amino acid mutations present in both subsets
   - Creates scatter plots of fitness effects with Pearson correlation
2. Filters to mutations meeting the count threshold in both subsets

### Step 17: Check Subset PCP Overlap

Checks whether parent-child pairs (PCPs) overlap between subsets within each grouping dimension:
1. Loads PCP files for each subset (host, geographic, temporal)
2. For each pair of subsets within a dimension, computes:
   - Parent-only overlap: fraction of parent sequences shared between subsets
   - Parent+child overlap: fraction of (parent, child) pairs shared between subsets
3. Reports counts and fractions, with the expectation that host and temporal subsets have small overlap, while geographic subsets may have substantial overlap

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

# Process geographic/temporal subset trees
snakemake --cores 8 results/HA/H3/north_america/mutation_counts.csv results/HA/H3/early/mutation_counts.csv

# Run subset fitness effect analysis
snakemake --cores 8 results/subset_aa_fitness_effects.csv
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
   - `parent_child_pairs.csv` - Detailed branch-level mutation information (same columns as global `parent_child_pairs.csv`)

3b. **Geographic Tree Outputs**: `results/{segment}/{subtype}/{geo}/`
   - `mutation_counts.csv` - Aggregated mutation counts for the specific geographic region (same columns as global `mutation_counts.csv`)
   - `parent_child_pairs.csv` - Detailed branch-level mutation information (same columns as global `parent_child_pairs.csv`)

3c. **Temporal Tree Outputs**: `results/{segment}/{subtype}/{temporal}/`
   - `mutation_counts.csv` - Aggregated mutation counts for the specific temporal group (same columns as global `mutation_counts.csv`)
   - `parent_child_pairs.csv` - Detailed branch-level mutation information (same columns as global `parent_child_pairs.csv`)

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
       - `host` — host group (`"human"`, `"avian"`, or `"all"`)
       - `evo_opp` — evolutionary opportunity
       - `rate` — per-site mutation rate
   - `genome_wide_rates.csv` - Mutation rates by type and class (synonymous, nonsynonymous, nonsense)
     - Columns:
       - `mut_type` — two-character mutation type (e.g. `"AC"` for A→C)
       - `mut_class` — mutation class (`"synonymous"`, `"nonsynonymous"`, or `"nonsense"`)
       - `host` — host group (`"human"`, `"avian"`, or `"all"`)
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
       - `host` — host group (`"human"`, `"avian"`, or `"all"`)
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
   - `nt_fitness_effects.csv` - Per-nucleotide-mutation fitness effects
     - Columns:
       - `host` — host group (`"human"`, `"avian"`, or `"all"`)
       - `subtype` — influenza subtype (e.g. `"H1"`, `"N2"`, or `"all"`)
       - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
       - `gene` — gene containing the site
       - `site` — nucleotide position in the reference sequence (1-indexed)
       - `wt_nt` — wildtype nucleotide
       - `mut_nt` — mutant nucleotide
       - `nt_mut` — nucleotide mutation string (e.g. `"A1C"`)
       - `mut_class` — mutation class (`"synonymous"`, `"nonsynonymous"`, `"nonsense"`, or `"noncoding"`)
       - `actual_count` — total observed occurrences of this nucleotide mutation
       - `expected_count` — total expected occurrences under the neutral model
       - `delta_fitness` — estimated fitness effect: log((actual_count + 0.5) / (expected_count + 0.5))
   - `sitewise_synonymous_fitness_effects.csv` - Per-site synonymous fitness effects
     - Columns:
       - `host` — host group (`"human"`, `"avian"`, or `"all"`)
       - `subtype` — influenza subtype (e.g. `"H1"`, `"N2"`, or `"all"`)
       - `segment` — genome segment (e.g. `"HA"`, `"PB2"`)
       - `gene` — gene containing the site
       - `site` — nucleotide position in the reference sequence (1-indexed)
       - `delta_fitness` — estimated fitness effect: mean of per-nucleotide-mutation `delta_fitness` values (from `nt_fitness_effects.csv`) across synonymous mutations at this site
   - `aa_fitness_effects.csv` - Per-amino-acid-mutation fitness effects
     - Columns:
       - `host` — host group (`"human"`, `"avian"`, or `"all"`)
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

9. **SHAPE-MaP Data**: `results/shapemap/`
   - `all_data.csv` - SHAPE-MaP RNA structure reactivities (Dadonaite et al. 2019) mapped to reference site coordinates for internal segments

10. **Subset Analysis Outputs**:
   - `subset_counts.csv` - Aggregated mutation counts from all subset trees (host, geographic, temporal)
     - Same columns as `counts.csv`, except `host` is replaced by:
       - `subset` — subset name (e.g. `"human"`, `"north_america"`, `"early"`)
       - `subset_type` — subset category (`"host"`, `"geographic"`, `"temporal"`)
   - `subset_aa_fitness_effects.csv` - Per-amino-acid-mutation fitness effects for each subset
     - Columns:
       - `subset` — subset name
       - `subset_type` — subset category
       - `subtype` — influenza subtype
       - `segment` — genome segment
       - `gene` — gene name
       - `codon_site` — codon site within the gene
       - `wt_aa` — wildtype amino acid
       - `mut_aa` — mutant amino acid
       - `aa_mut` — amino acid mutation string
       - `mut_class` — mutation class
       - `actual_count` — total observed mutations
       - `expected_count` — total expected mutations under the neutral model
       - `delta_fitness` — estimated fitness effect: log((actual_count + 0.5) / (expected_count + 0.5))

11. **Processed DMS Data**: `results/dms_data/`
   - `Yu_HA/processed_dms_data.csv` - HA DMS phenotypes (Yu et al. 2025) with tree reference site numbering
   - `Bloom_NP/processed_dms_data.csv` - NP DMS preferences (Bloom et al. 2014) with log-ratio fitness effects
   - `.process_dms_data_soh_pb2.done` - Sentinel file marking completion of Soh et al. PB2 alignment QC
   - `Wang_NA/processed_dms_data.csv` - NA DMS fitness effects (Wang et al. 2023) with tree reference site numbering
   - `Li_PB1/processed_dms_data.csv` - PB1 DMS fitness effects (Li et al. 2023) with tree reference site numbering
   - `Hom_M1/processed_dms_data.csv` - M1 DMS fitness effects (Hom et al. 2019) with tree reference site numbering
   - `Teo_NEP/processed_dms_data.csv` - NEP DMS fitness effects (Teo et al. 2024) with tree reference site numbering
   - `Chen_PA/processed_dms_data.csv` - PA DMS fitness effects (Chen et al. 2024) for first 240 AA with tree reference site numbering

### Output Structure Example

```
results/
├── HA/
│   ├── H1/
│   │   ├── coding_sites.csv
│   │   ├── mutation_counts.csv
│   │   ├── parent_child_pairs.csv
│   │   ├── human/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   ├── avian/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   ├── north_america/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   ├── europe/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   ├── asia/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   ├── early/
│   │   │   ├── mutation_counts.csv
│   │   │   └── parent_child_pairs.csv
│   │   └── late/
│   │       ├── mutation_counts.csv
│   │       └── parent_child_pairs.csv
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
├── nt_fitness_effects.csv
├── sitewise_synonymous_fitness_effects.csv
├── aa_fitness_effects.csv
├── shapemap/
│   └── all_data.csv
├── dms_data/
│   ├── Yu_HA/
│   │   └── processed_dms_data.csv
│   ├── Bloom_NP/
│   │   └── processed_dms_data.csv
│   ├── Wang_NA/
│   │   └── processed_dms_data.csv
│   ├── Li_PB1/
│   │   └── processed_dms_data.csv
│   ├── Hom_M1/
│   │   └── processed_dms_data.csv
│   ├── Teo_NEP/
│   │   └── processed_dms_data.csv
│   └── Chen_PA/
│       └── processed_dms_data.csv
├── .process_dms_data_soh_pb2.done
├── .summarize_filter_logs.done
├── subset_counts.csv
├── subset_aa_fitness_effects.csv
├── .analyze_subset_fitness_effects.done
├── .check_subset_pcp_overlap.done
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

## Interactive Dashboards

Two interactive fitness effect dashboards are hosted via GitHub Pages under `docs/` and exported using [marimo](https://marimo.io/) in WASM mode:

- **AA dashboard** (`docs/aa/index.html`): per-amino-acid-mutation fitness effects (`notebooks/aa_fitness_heatmap_dashboard.py`)
- **NT dashboard** (`docs/nt/index.html`): per-nucleotide-mutation fitness effects (`notebooks/nt_fitness_heatmap_dashboard.py`)

To regenerate the HTML exports after modifying a dashboard notebook, run the corresponding Snakemake rule:

```bash
# Regenerate the NT dashboard
snakemake --cores 1 docs/nt/index.html --forcerun export_nt_dashboard

# Regenerate the AA dashboard
snakemake --cores 1 docs/aa/index.html --forcerun export_dashboard
```

Each rule exports the notebook as a self-contained WASM app, patches the theme to dark, and copies the required result files into the `docs/{aa,nt}/results/` directory so they are accessible to the app at runtime.

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
conda activate flu-mut-rates
```

The pipeline will automatically execute Jupyter notebooks as part of the workflow using `jupyter nbconvert`.

