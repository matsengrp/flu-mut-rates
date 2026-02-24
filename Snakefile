# Load configuration
configfile: "config.yaml"

# Helper function to get subtypes for each segment
def get_subtypes_for_segment(segment):
    if segment == "HA":
        return config["ha_subtypes"]
    elif segment == "NA":
        return config["na_subtypes"]
    else:
        return ["all"]

def get_all_mutation_count_files():
    """Generate list of all mutation_counts.csv files (global + host-specific)"""
    files = []

    # Global mutation counts (no host subdirectory)
    for segment in config["segments"]:
        for subtype in get_subtypes_for_segment(segment):
            files.append(
                f"{config['output_dir']}/{segment}/{subtype}/mutation_counts.csv"
            )

    # Host-specific mutation counts (human and avian only, excluding swine)
    for segment in config["segments"]:
        for subtype in get_subtypes_for_segment(segment):
            for host in config["host_groups"]:
                files.append(
                    f"{config['output_dir']}/{segment}/{subtype}/{host}/mutation_counts.csv"
                )

    return files

# Generate all segment-subtype-host combinations for host-specific trees
segment_subtype_host_combinations = []
for segment in config["segments"]:
    for subtype in get_subtypes_for_segment(segment):
        for host in config["host_groups"]:
            segment_subtype_host_combinations.append((segment, subtype, host))

# Generate all segment-subtype combinations for global trees
segment_subtype_combinations = []
for segment in config["segments"]:
    for subtype in get_subtypes_for_segment(segment):
        segment_subtype_combinations.append((segment, subtype))

# Final output files
final_outputs = []

# Add host-specific mutation counts (no PCPs)
for segment, subtype, host in segment_subtype_host_combinations:
    final_outputs.append(
        f"{config['output_dir']}/{segment}/{subtype}/{host}/mutation_counts.csv"
    )

# Add global mutation counts and PCPs
for segment, subtype in segment_subtype_combinations:
    final_outputs.extend([
        f"{config['output_dir']}/{segment}/{subtype}/mutation_counts.csv",
        f"{config['output_dir']}/{segment}/{subtype}/parent_child_pairs.csv"
    ])

# Add aligned proteins outputs (only for HA and NA segments)
final_outputs.extend(expand(
    "{output_dir}/aligned_proteins/{segment}",
    output_dir=config["output_dir"],
    segment=["HA", "NA"]
))

# Add compute_rates outputs to final targets
final_outputs.extend([
    f"{config['output_dir']}/counts.csv",
    f"{config['output_dir']}/genome_wide_rates.csv",
    f"{config['output_dir']}/motif_level_genome_wide_rates.csv",
    f"{config['output_dir']}/evo_opp_thresholds.csv",
    f"{config['output_dir']}/site_specific_mutation_rates.csv"
])

# Add genome-wide rates analysis notebook to final targets
final_outputs.extend([
    f"{config['output_dir']}/.analyze_genome_wide_rates.done"
])

# Add neutral model outputs to final targets
final_outputs.extend([
    f"{config['output_dir']}/neutral_model/base/model_performance.csv",
    f"{config['output_dir']}/neutral_model/local_context/model_performance.csv",
    f"{config['output_dir']}/neutral_model/local_context+global_context/model_performance.csv",
    f"{config['output_dir']}/expected_rates.csv"
])

# Main rule to define target outputs
rule all:
    input:
        final_outputs

# Create coding sites file
rule make_coding_sites:
    input:
        ref_fasta=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gff_file=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gff"
    output:
        coding_sites="{output_dir}/{segment}/{subtype}/coding_sites.csv"
    log:
        "logs/{output_dir}/{segment}_{subtype}_coding_sites.log"
    shell:
        """
        python scripts/make_coding_sites.py \
            --ref_fasta {input.ref_fasta} \
            --gff_file {input.gff_file} \
            --subtype {wildcards.subtype} \
            --segment {wildcards.segment} \
            --output_directory {wildcards.output_dir}/{wildcards.segment}/{wildcards.subtype} &> {log}
        """

# Count mutations along host-specific trees
rule count_mutations_host_trees:
    input:
        tree_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/host_specific_trees/{wildcards.host}_tree.pb.gz",
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gtf_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gtf"
    output:
        all_counts_path="{output_dir}/{segment}/{subtype}/{host}/mutation_counts.csv"
    log:
        "logs/{output_dir}/{segment}_{subtype}_{host}_mutation_counts.log"
    shell:
        """
        python scripts/make_count_dfs.py \
            --tree_path {input.tree_path} \
            --coding_site_path {input.coding_site_path} \
            --fasta_path {input.fasta_path} \
            --gtf_path {input.gtf_path} \
            --all_counts_path {output.all_counts_path} &> {log}
        """

# Count mutations along global (non-host-specific) trees
rule count_mutations:
    input:
        tree_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/final_tree.pb.gz",
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gtf_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gtf"
    output:
        all_counts_path="{output_dir}/{segment}/{subtype}/mutation_counts.csv",
        all_pcps_path="{output_dir}/{segment}/{subtype}/parent_child_pairs.csv"
    log:
        "logs/{output_dir}/{segment}_{subtype}_mutation_counts.log"
    shell:
        """
        python scripts/make_count_dfs.py \
            --tree_path {input.tree_path} \
            --coding_site_path {input.coding_site_path} \
            --fasta_path {input.fasta_path} \
            --gtf_path {input.gtf_path} \
            --all_counts_path {output.all_counts_path} \
            --all_pcps_path {output.all_pcps_path} &> {log}
        """

# Compute genome-wide and site-specific mutation rates from all mutation count data
rule compute_rates:
    input:
        notebook="notebooks/compute_rates.ipynb",
        mutation_counts=get_all_mutation_count_files()
    output:
        counts="{output_dir}/counts.csv",
        genome_wide_rates="{output_dir}/genome_wide_rates.csv",
        motif_rates="{output_dir}/motif_level_genome_wide_rates.csv",
        evo_opp_thresholds="{output_dir}/evo_opp_thresholds.csv",
        site_specific_rates="{output_dir}/site_specific_mutation_rates.csv"
    log:
        "logs/{output_dir}/compute_rates.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            compute_rates.ipynb &> ../{log}
        """

# Analyze genome-wide mutation rates and create visualizations
rule analyze_genome_wide_rates:
    input:
        notebook="notebooks/analyze_genome_wide_rates.ipynb",
        genome_wide_rates="{output_dir}/genome_wide_rates.csv"
    output:
        touch("{output_dir}/.analyze_genome_wide_rates.done")
    log:
        "logs/{output_dir}/analyze_genome_wide_rates.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            analyze_genome_wide_rates.ipynb &> ../{log}
        """

# Fit neutral linear models to synonymous mutation rates
rule fit_neutral_models:
    input:
        script="scripts/rates_model.py",
        site_specific_rates="{output_dir}/site_specific_mutation_rates.csv"
    output:
        # Base model outputs (no factors)
        base_expected="{output_dir}/neutral_model/base/expected_rates_by_predictor.csv",
        base_performance="{output_dir}/neutral_model/base/model_performance.csv",
        # Local context model outputs (motif only)
        local_expected="{output_dir}/neutral_model/local_context/expected_rates_by_predictor.csv",
        local_performance="{output_dir}/neutral_model/local_context/model_performance.csv",
        # Local + global context model outputs (motif + segment)
        full_expected="{output_dir}/neutral_model/local_context+global_context/expected_rates_by_predictor.csv",
        full_performance="{output_dir}/neutral_model/local_context+global_context/model_performance.csv"
    log:
        "logs/{output_dir}/fit_neutral_models.log"
    shell:
        """
        python {input.script} &> {log}
        """

# Augment neutral model expected rates with CG and GC mutation types
rule augment_expected_rates:
    input:
        script="scripts/augment_expected_rates.py",
        segment_wide_rates="{output_dir}/segment_wide_rates.csv",
        full_expected="{output_dir}/neutral_model/local_context+global_context/expected_rates_by_predictor.csv"
    output:
        expected_rates="{output_dir}/expected_rates.csv"
    log:
        "logs/{output_dir}/augment_expected_rates.log"
    shell:
        """
        python {input.script} \
            --segment_wide_rates {input.segment_wide_rates} \
            --input_file {input.full_expected} \
            --output_file {output.expected_rates} &> {log}
        """

# Align protein sequences across subtypes (only for HA and NA)
rule align_proteins:
    input:
        # Ensure all coding sites files are created first for this segment
        coding_sites=lambda wildcards: expand("{output_dir}/{segment}/{subtype}/coding_sites.csv",
                          output_dir=config["output_dir"],
                          segment=wildcards.segment,
                          subtype=get_subtypes_for_segment(wildcards.segment))
    output:
        # Output directory marker (you could also specify specific aligned files)
        aligned_dir=directory("{output_dir}/aligned_proteins/{segment}")
    log:
        "logs/{output_dir}/align_proteins_{segment}.log"
    wildcard_constraints:
        segment="HA|NA"
    params:
        subtypes=lambda wildcards: get_subtypes_for_segment(wildcards.segment)
    shell:
        """
        python scripts/align_proteins.py \
            --output_dir {config[output_dir]} \
            --segment {wildcards.segment} \
            --subtypes {params.subtypes} \
            --muscle_path muscle &> {log}
        """