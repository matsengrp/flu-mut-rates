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

# Segments used in SHAPE-MaP analysis (HA and NA excluded per the notebook)
shapemap_segments = [s for s in config["segments"] if s not in ["HA", "NA"]]

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

# Add SHAPE-MaP processed data to final targets
final_outputs.append(f"{config['output_dir']}/shapemap/all_data.csv")

# Add site-specific rates analysis notebook to final targets
final_outputs.extend([
    f"{config['output_dir']}/.analyze_site_specific_rates.done"
])

# Add compute_fitness_effects outputs to final targets
final_outputs.extend([
    f"{config['output_dir']}/actual_expected.csv",
    f"{config['output_dir']}/sitewise_synonymous_fitness_effects.csv",
    f"{config['output_dir']}/aa_fitness_effects.csv"
])

# Add process_dms_data outputs to final targets
final_outputs.extend([
    f"{config['output_dir']}/dms_data/Yu_HA/processed_dms_data.csv",
    f"{config['output_dir']}/dms_data/Bloom_NP/processed_dms_data.csv",
    f"{config['output_dir']}/.process_dms_data_soh_pb2.done",
    f"{config['output_dir']}/.process_dms_data_wang_na.done",
])

# Add analyze_fitness_effects to final targets
final_outputs.append(f"{config['output_dir']}/.analyze_fitness_effects.done")

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
        "logs/{segment}/{subtype}/coding_sites.log"
    params:
        ignore_genes=lambda wildcards: " ".join(
            config.get("ignore_genes", {}).get(wildcards.segment, [])
        )
    shell:
        """
        python scripts/make_coding_sites.py \
            --ref_fasta {input.ref_fasta} \
            --gff_file {input.gff_file} \
            --subtype {wildcards.subtype} \
            --segment {wildcards.segment} \
            --output_directory {wildcards.output_dir}/{wildcards.segment}/{wildcards.subtype} \
            --ignore_genes {params.ignore_genes} &> {log}
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
        "logs/{segment}/{subtype}/{host}/mutation_counts.log"
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
        "logs/{segment}/{subtype}/mutation_counts.log"
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
        "logs/compute_rates.log"
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
        "logs/analyze_genome_wide_rates.log"
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

# Analyze site-specific mutation rates and create visualizations
rule analyze_site_specific_rates:
    input:
        notebook="notebooks/analyze_site_specific_rates.ipynb",
        site_specific_rates="{output_dir}/site_specific_mutation_rates.csv",
        segment_wide_rates="{output_dir}/segment_wide_rates.csv",
        motif_rates="{output_dir}/motif_level_genome_wide_rates.csv",
        local_context_performance="{output_dir}/neutral_model/local_context/model_performance.csv",
        full_context_performance="{output_dir}/neutral_model/local_context+global_context/model_performance.csv",
        shapemap_data="{output_dir}/shapemap/all_data.csv",
        motif_medians="data/haddox_2025/motif_medians.csv"
    output:
        touch("{output_dir}/.analyze_site_specific_rates.done")
    log:
        "logs/analyze_site_specific_rates.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            analyze_site_specific_rates.ipynb &> ../{log}
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
        "logs/fit_neutral_models.log"
    shell:
        """
        python {input.script} &> {log}
        """

# Augment neutral model expected rates with CG and GC mutation types
rule augment_expected_rates:
    input:
        script="scripts/augment_expected_rates.py",
        segment_wide_rates="{output_dir}/segment_wide_rates.csv",
        motif_level_rates="{output_dir}/motif_level_genome_wide_rates.csv",
        full_expected="{output_dir}/neutral_model/local_context+global_context/expected_rates_by_predictor.csv"
    output:
        expected_rates="{output_dir}/expected_rates.csv"
    log:
        "logs/augment_expected_rates.log"
    shell:
        """
        python {input.script} \
            --segment_wide_rates {input.segment_wide_rates} \
            --motif_level_rates {input.motif_level_rates} \
            --input_file {input.full_expected} \
            --output_file {output.expected_rates} &> {log}
        """

# Compute fitness effects of synonymous and amino-acid mutations
rule compute_fitness_effects:
    input:
        notebook="notebooks/compute_fitness_effects.ipynb",
        counts="{output_dir}/counts.csv",
        expected_rates="{output_dir}/expected_rates.csv"
    output:
        actual_expected="{output_dir}/actual_expected.csv",
        syn_fitness="{output_dir}/sitewise_synonymous_fitness_effects.csv",
        aa_fitness="{output_dir}/aa_fitness_effects.csv"
    log:
        "logs/compute_fitness_effects.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            compute_fitness_effects.ipynb &> ../{log}
        """

# Process HA DMS data from Yu et al.
rule process_dms_data_yu_ha:
    input:
        notebook="notebooks/process_dms_data_yu_ha.ipynb",
        ha_phenotypes="data/dms_data/Yu_HA/Phenotypes.csv",
        ha_numbering="data/dms_data/Yu_HA/site_numbering_map.csv"
    output:
        ha_dms="{output_dir}/dms_data/Yu_HA/processed_dms_data.csv"
    log:
        "logs/process_dms_data_yu_ha.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_yu_ha.ipynb &> ../{log}
        """

# Process NP DMS data from Bloom et al.
rule process_dms_data_bloom_np:
    input:
        notebook="notebooks/process_dms_data_bloom_np.ipynb",
        np_data="data/dms_data/Bloom_NP/Supplementary_file_1.xls"
    output:
        np_dms="{output_dir}/dms_data/Bloom_NP/processed_dms_data.csv"
    log:
        "logs/process_dms_data_bloom_np.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_bloom_np.ipynb &> ../{log}
        """

# Process PB2 DMS data from Soh et al. (alignment QC)
rule process_dms_data_soh_pb2:
    input:
        notebook="notebooks/process_dms_data_soh_pb2.ipynb",
        pb2_data="data/dms_data/Soh_PB2/elife-45079-fig2-data1-v1.csv"
    output:
        touch("{output_dir}/.process_dms_data_soh_pb2.done")
    log:
        "logs/process_dms_data_soh_pb2.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_soh_pb2.ipynb &> ../{log}
        """

# Process NA DMS data from Wang et al. (sequence comparison)
rule process_dms_data_wang_na:
    input:
        notebook="notebooks/process_dms_data_wang_na.ipynb",
        na_data="data/dms_data/Wang_NA/msystems.00670-23-s0006.xlsx"
    output:
        touch("{output_dir}/.process_dms_data_wang_na.done")
    log:
        "logs/process_dms_data_wang_na.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_wang_na.ipynb &> ../{log}
        """

# Analyze fitness effects and compare to DMS data
rule analyze_fitness_effects:
    input:
        notebook="notebooks/analyze_fitness_effects.ipynb",
        aa_fitness="{output_dir}/aa_fitness_effects.csv",
        syn_fitness="{output_dir}/sitewise_synonymous_fitness_effects.csv",
        ha_dms="{output_dir}/dms_data/Yu_HA/processed_dms_data.csv",
        np_dms="{output_dir}/dms_data/Bloom_NP/processed_dms_data.csv",
        pb2_dms="data/dms_data/Soh_PB2/elife-45079-fig2-data1-v1.csv"
    output:
        touch("{output_dir}/.analyze_fitness_effects.done")
    log:
        "logs/analyze_fitness_effects.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            analyze_fitness_effects.ipynb &> ../{log}
        """

# Process SHAPE-MaP reactivity data and align to reference sequences
rule process_shapemap_data:
    input:
        excel="data/dadonaite_2019/41564_2019_513_MOESM3_ESM.xlsx",
        ref_fastas=expand(
            "{data_dir}/{segment}/all/curated_reference.fasta",
            data_dir=config["data_dir"],
            segment=shapemap_segments
        )
    output:
        "{output_dir}/shapemap/all_data.csv"
    params:
        data_dir=config["data_dir"]
    log:
        "logs/process_shapemap_data.log"
    shell:
        """
        data_dir=$(realpath {params.data_dir}) && \
        cd notebooks && \
        papermill \
            process_shapemap_data.ipynb \
            process_shapemap_data.ipynb \
            -p data_dir $data_dir &> ../{log}
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
        "logs/{segment}/align_proteins.log"
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