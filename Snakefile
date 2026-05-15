import json

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
    """Generate list of mutation_counts.csv files for the global tree and the
    host-stratified output (one file per segment/subtype, with a host column)."""
    files = []

    # Global mutation counts (no subdirectory)
    for segment in config["segments"]:
        for subtype in get_subtypes_for_segment(segment):
            files.append(
                f"{config['output_dir']}/{segment}/{subtype}/mutation_counts.csv"
            )

    # Host-stratified mutation counts (single file per tree with a host column)
    for segment in config["segments"]:
        for subtype in get_subtypes_for_segment(segment):
            files.append(
                f"{config['output_dir']}/{segment}/{subtype}/host_stratified/mutation_counts.csv"
            )

    return files

def get_subset_mutation_count_files():
    """Generate list of subset mutation_counts.csv files for downstream subset rates.
    The host dimension is supplied by the single host_stratified file (with a host
    column); geographic and split-half dimensions come from per-subset files
    (split-half files are written by the global count_mutations rule in the same
    pass)."""
    files = []

    for segment in config["segments"]:
        for subtype in get_subtypes_for_segment(segment):
            files.append(
                f"{config['output_dir']}/{segment}/{subtype}/host_stratified/mutation_counts.csv"
            )
            for geo in config["geographic_groups"]:
                files.append(
                    f"{config['output_dir']}/{segment}/{subtype}/{geo}/mutation_counts.csv"
                )
            for split in config["split_half_groups"]:
                files.append(
                    f"{config['output_dir']}/{segment}/{subtype}/{split}/mutation_counts.csv"
                )

    return files

# Generate all segment-subtype-geographic combinations for geographic trees
segment_subtype_geo_combinations = []
for segment in config["segments"]:
    for subtype in get_subtypes_for_segment(segment):
        for geo in config["geographic_groups"]:
            segment_subtype_geo_combinations.append((segment, subtype, geo))

# Generate all segment-subtype combinations for global trees
segment_subtype_combinations = []
for segment in config["segments"]:
    for subtype in get_subtypes_for_segment(segment):
        segment_subtype_combinations.append((segment, subtype))

# Segments used in SHAPE-MaP analysis (HA and NA excluded per the notebook)
shapemap_segments = [s for s in config["segments"] if s not in ["HA", "NA"]]

# Final output files
final_outputs = []

# Add geographic mutation counts and PCPs
for segment, subtype, geo in segment_subtype_geo_combinations:
    final_outputs.extend([
        f"{config['output_dir']}/{segment}/{subtype}/{geo}/mutation_counts.csv",
        f"{config['output_dir']}/{segment}/{subtype}/{geo}/parent_child_pairs.csv"
    ])

# Add global mutation counts and PCPs
for segment, subtype in segment_subtype_combinations:
    final_outputs.extend([
        f"{config['output_dir']}/{segment}/{subtype}/mutation_counts.csv",
        f"{config['output_dir']}/{segment}/{subtype}/parent_child_pairs.csv"
    ])

# Add host-stratified mutation counts and PCPs (global tree + host TSV)
for segment, subtype in segment_subtype_combinations:
    final_outputs.extend([
        f"{config['output_dir']}/{segment}/{subtype}/host_stratified/mutation_counts.csv",
        f"{config['output_dir']}/{segment}/{subtype}/host_stratified/parent_child_pairs.csv"
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
    f"{config['output_dir']}/site_specific_substitution_rates.csv"
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
    f"{config['output_dir']}/sitewise_synonymous_fitness_effects.csv",
    f"{config['output_dir']}/aa_fitness_effects.csv"
])

# Add process_dms_data outputs to final targets
final_outputs.extend([
    f"{config['output_dir']}/dms_data/Yu_HA/processed_dms_data.csv",
    f"{config['output_dir']}/dms_data/Bloom_NP/processed_dms_data.csv",
    f"{config['output_dir']}/.process_dms_data_soh_pb2.done",
    f"{config['output_dir']}/dms_data/Wang_NA/processed_dms_data.csv",
    f"{config['output_dir']}/dms_data/Li_PB1/processed_dms_data.csv",
    f"{config['output_dir']}/dms_data/Hom_M1/processed_dms_data.csv",
    f"{config['output_dir']}/dms_data/Teo_NEP/processed_dms_data.csv",
    f"{config['output_dir']}/dms_data/Chen_PA/processed_dms_data.csv",
])

# Add analyze_fitness_effects to final targets
final_outputs.append(f"{config['output_dir']}/.analyze_fitness_effects.done")

# Add summarize_filter_logs notebook to final targets
final_outputs.append(f"{config['output_dir']}/.summarize_filter_logs.done")

# Add subset analysis outputs to final targets
final_outputs.extend([
    f"{config['output_dir']}/subset_counts.csv",
    f"{config['output_dir']}/subset_aa_fitness_effects.csv",
    f"{config['output_dir']}/.analyze_subset_fitness_effects.done",
])

# Add compose_figures notebook to final targets
final_outputs.append(f"{config['output_dir']}/.compose_figures.done")

# Add dashboard exports to final targets
final_outputs.append("docs/index.html")
final_outputs.append("docs/aa/index.html")
final_outputs.append("docs/nt/index.html")

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
        coding_sites="{output_dir}/{segment}/{subtype}/coding_sites.csv",
        protein_sequences=directory("{output_dir}/{segment}/{subtype}/protein_sequences")
    log:
        "{output_dir}/logs/{segment}/{subtype}/coding_sites.log"
    params:
        ignore_genes=lambda wildcards: " ".join(
            config.get("ignore_genes", {}).get(wildcards.segment, [])
        ),
        additional_orfs=lambda wildcards: json.dumps(
            config.get("additional_orfs", {}).get(wildcards.segment, [])
        )
    shell:
        """
        python scripts/make_coding_sites.py \
            --ref_fasta {input.ref_fasta} \
            --gff_file {input.gff_file} \
            --subtype {wildcards.subtype} \
            --segment {wildcards.segment} \
            --output_directory {wildcards.output_dir}/{wildcards.segment}/{wildcards.subtype} \
            --ignore_genes {params.ignore_genes} \
            --additional_orfs '{params.additional_orfs}' &> {log}
        """

# Count mutations along geographic trees
rule count_mutations_geographic_trees:
    input:
        tree_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/geographic_trees/{wildcards.geo}_tree.pb.gz",
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gtf_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gtf"
    output:
        all_counts_path="{output_dir}/{segment}/{subtype}/{geo}/mutation_counts.csv",
        all_pcps_path="{output_dir}/{segment}/{subtype}/{geo}/parent_child_pairs.csv"
    log:
        "{output_dir}/logs/{segment}/{subtype}/{geo}/mutation_counts.log"
    wildcard_constraints:
        geo="|".join(config["geographic_groups"])
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

# Count mutations along the global tree, stratified by host using ancestral-state TSV.
# Only branches whose parent and child share the same unambiguous host are counted,
# and the output is restricted to the host labels listed in config["host_groups"].
rule count_mutations_host_stratified:
    input:
        tree_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/final_tree.pb.gz",
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gtf_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gtf",
        host_tsv=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/host_ancestral/combined_ancestral_states.tab"
    output:
        all_counts_path="{output_dir}/{segment}/{subtype}/host_stratified/mutation_counts.csv",
        all_pcps_path="{output_dir}/{segment}/{subtype}/host_stratified/parent_child_pairs.csv"
    log:
        "{output_dir}/logs/{segment}/{subtype}/host_stratified/mutation_counts.log"
    params:
        host_groups=" ".join(config["host_groups"])
    shell:
        """
        python scripts/make_count_dfs.py \
            --tree_path {input.tree_path} \
            --coding_site_path {input.coding_site_path} \
            --fasta_path {input.fasta_path} \
            --gtf_path {input.gtf_path} \
            --host_tsv {input.host_tsv} \
            --host_groups {params.host_groups} \
            --all_counts_path {output.all_counts_path} \
            --all_pcps_path {output.all_pcps_path} &> {log}
        """

# Count mutations along global (non-host-specific) trees. The same traversal also
# emits two complementary random-split-half mutation_counts files (split_a, split_b)
# for use as a fitness-effect noise floor; the seed is fixed in config.yaml so
# results are reproducible.
rule count_mutations:
    input:
        tree_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/final_tree.pb.gz",
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gtf_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gtf"
    output:
        all_counts_path="{output_dir}/{segment}/{subtype}/mutation_counts.csv",
        all_pcps_path="{output_dir}/{segment}/{subtype}/parent_child_pairs.csv",
        all_counts_path_split_a="{output_dir}/{segment}/{subtype}/split_a/mutation_counts.csv",
        all_counts_path_split_b="{output_dir}/{segment}/{subtype}/split_b/mutation_counts.csv"
    log:
        "{output_dir}/logs/{segment}/{subtype}/mutation_counts.log"
    params:
        partition_seed=config["split_half_seed"]
    shell:
        """
        python scripts/make_count_dfs.py \
            --tree_path {input.tree_path} \
            --coding_site_path {input.coding_site_path} \
            --fasta_path {input.fasta_path} \
            --gtf_path {input.gtf_path} \
            --all_counts_path {output.all_counts_path} \
            --all_pcps_path {output.all_pcps_path} \
            --partition_seed {params.partition_seed} \
            --all_counts_path_split_a {output.all_counts_path_split_a} \
            --all_counts_path_split_b {output.all_counts_path_split_b} &> {log}
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
        site_specific_rates="{output_dir}/site_specific_substitution_rates.csv",
        segment_wide_rates="{output_dir}/segment_wide_rates.csv"
    log:
        "{output_dir}/logs/compute_rates.log"
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
        "{output_dir}/logs/analyze_genome_wide_rates.log"
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
        site_specific_rates="{output_dir}/site_specific_substitution_rates.csv",
        segment_wide_rates="{output_dir}/segment_wide_rates.csv",
        motif_rates="{output_dir}/motif_level_genome_wide_rates.csv",
        local_context_performance="{output_dir}/neutral_model/local_context/model_performance.csv",
        full_context_performance="{output_dir}/neutral_model/local_context+global_context/model_performance.csv",
        shapemap_data="{output_dir}/shapemap/all_data.csv",
        motif_medians="data/haddox_2025/motif_medians.csv"
    output:
        touch("{output_dir}/.analyze_site_specific_rates.done")
    log:
        "{output_dir}/logs/analyze_site_specific_rates.log"
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
        site_specific_rates="{output_dir}/site_specific_substitution_rates.csv"
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
        "{output_dir}/logs/fit_neutral_models.log"
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
        "{output_dir}/logs/augment_expected_rates.log"
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
        syn_fitness="{output_dir}/sitewise_synonymous_fitness_effects.csv",
        aa_fitness="{output_dir}/aa_fitness_effects.csv",
        nt_fitness="{output_dir}/nt_fitness_effects.csv"
    log:
        "{output_dir}/logs/compute_fitness_effects.log"
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
        "{output_dir}/logs/process_dms_data_yu_ha.log"
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
        "{output_dir}/logs/process_dms_data_bloom_np.log"
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
        "{output_dir}/logs/process_dms_data_soh_pb2.log"
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
        na_dms="{output_dir}/dms_data/Wang_NA/processed_dms_data.csv"
    log:
        "{output_dir}/logs/process_dms_data_wang_na.log"
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

# Process PB1 DMS data from Li et al.
rule process_dms_data_li_pb1:
    input:
        notebook="notebooks/process_dms_data_li_pb1.ipynb",
        pb1_data="data/dms_data/Li_PB1/jvi.01329-23-s0008.csv"
    output:
        pb1_dms="{output_dir}/dms_data/Li_PB1/processed_dms_data.csv"
    log:
        "{output_dir}/logs/process_dms_data_li_pb1.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_li_pb1.ipynb &> ../{log}
        """

# Process M1 DMS data from Hom et al.
rule process_dms_data_hom_m1:
    input:
        notebook="notebooks/process_dms_data_hom_m1.ipynb",
        prefs="data/dms_data/Hom_M1/summary_avgprefs.csv",
        fasta="data/dms_data/Hom_M1/PR8-M1.fasta"
    output:
        m1_dms="{output_dir}/dms_data/Hom_M1/processed_dms_data.csv"
    log:
        "{output_dir}/logs/process_dms_data_hom_m1.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_hom_m1.ipynb &> ../{log}
        """

# Process NEP DMS data from Teo et al.
rule process_dms_data_teo_nep:
    input:
        notebook="notebooks/process_dms_data_teo_nep.ipynb",
        nep_data="data/dms_data/Teo_NEP/mmc2.xlsx"
    output:
        nep_dms="{output_dir}/dms_data/Teo_NEP/processed_dms_data.csv"
    log:
        "{output_dir}/logs/process_dms_data_teo_nep.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_teo_nep.ipynb &> ../{log}
        """

# Process PA DMS data from Chen et al.
rule process_dms_data_chen_pa:
    input:
        notebook="notebooks/process_dms_data_chen_pa.ipynb",
        pa_data="data/dms_data/Chen_PA/fitness calculation.xlsx"
    output:
        pa_dms="{output_dir}/dms_data/Chen_PA/processed_dms_data.csv"
    log:
        "{output_dir}/logs/process_dms_data_chen_pa.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_dms_data_chen_pa.ipynb &> ../{log}
        """


# Analyze fitness effects and compare to DMS data
rule analyze_fitness_effects:
    input:
        notebook="notebooks/analyze_fitness_effects.ipynb",
        aa_fitness="{output_dir}/aa_fitness_effects.csv",
        syn_fitness="{output_dir}/sitewise_synonymous_fitness_effects.csv",
        ha_dms="{output_dir}/dms_data/Yu_HA/processed_dms_data.csv",
        np_dms="{output_dir}/dms_data/Bloom_NP/processed_dms_data.csv",
        pb2_dms="data/dms_data/Soh_PB2/elife-45079-fig2-data1-v1.csv",
        na_dms="{output_dir}/dms_data/Wang_NA/processed_dms_data.csv",
        pb1_dms="{output_dir}/dms_data/Li_PB1/processed_dms_data.csv",
        m1_dms="{output_dir}/dms_data/Hom_M1/processed_dms_data.csv",
        nep_dms="{output_dir}/dms_data/Teo_NEP/processed_dms_data.csv",
        pa_dms="{output_dir}/dms_data/Chen_PA/processed_dms_data.csv"
    output:
        touch("{output_dir}/.analyze_fitness_effects.done")
    log:
        "{output_dir}/logs/analyze_fitness_effects.log"
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
    log:
        "{output_dir}/logs/process_shapemap_data.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            process_shapemap_data.ipynb &> ../{log}
        """

# Align protein sequences across subtypes (only for HA and NA)
rule align_proteins:
    input:
        # Ensure all coding sites files are created first for this segment
        coding_sites=lambda wildcards: expand("{output_dir}/{segment}/{subtype}/coding_sites.csv",
                          output_dir=config["output_dir"],
                          segment=wildcards.segment,
                          subtype=get_subtypes_for_segment(wildcards.segment)),
        protein_sequences=lambda wildcards: expand("{output_dir}/{segment}/{subtype}/protein_sequences",
                          output_dir=config["output_dir"],
                          segment=wildcards.segment,
                          subtype=get_subtypes_for_segment(wildcards.segment))
    output:
        # Output directory marker (you could also specify specific aligned files)
        aligned_dir=directory("{output_dir}/aligned_proteins/{segment}")
    log:
        "{output_dir}/logs/{segment}/align_proteins.log"
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

# Pre-compute reference amino acid sequences for all gene/subtype combinations
rule precompute_reference_aa:
    input:
        script="scripts/precompute_reference_aa.py"
    output:
        "{output_dir}/reference_aa.json"
    log:
        "{output_dir}/logs/precompute_reference_aa.log"
    shell:
        """
        python {input.script} &> {log}
        """

# Pre-compute reference nucleotide sequences for all segment/subtype combinations
rule precompute_reference_nt:
    input:
        script="scripts/precompute_reference_nt.py"
    output:
        "{output_dir}/reference_nt.json"
    log:
        "{output_dir}/logs/precompute_reference_nt.log"
    shell:
        """
        python {input.script} &> {log}
        """

# Landing page linking to both dashboards
rule create_landing_page:
    output:
        "docs/index.html"
    shell:
        """
        cat > docs/index.html << 'EOF'
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Flu Fitness Dashboards</title>
  <style>
    body {{ font-family: sans-serif; background: #121212; color: #e0e0e0;
           display: flex; flex-direction: column; align-items: center;
           justify-content: center; height: 100vh; margin: 0; }}
    h1 {{ margin-bottom: 2rem; }}
    .links {{ display: flex; gap: 2rem; }}
    a {{ display: block; padding: 1rem 2rem; background: #1e1e1e;
         border: 1px solid #BB86FC; border-radius: 8px; color: #BB86FC;
         text-decoration: none; font-size: 1.2rem; text-align: center; }}
    a:hover {{ background: #2a2a2a; }}
  </style>
</head>
<body>
  <h1>Flu Fitness Dashboards</h1>
  <div class="links">
    <a href="aa/">Amino acid mutations</a>
    <a href="nt/">Nucleotide mutations</a>
  </div>
  <p style="margin-top: 2rem; font-size: 0.9rem; color: #aaa; max-width: 500px; text-align: center; line-height: 1.5;">
    Fitness effects were estimated in Haddox et al. Interactive plots were made
    using code adapted from the <a href="https://jbloomlab.github.io/polyclonal/" style="display: inline; color: #BB86FC; border: none; padding: 0; background: none; font-size: inherit;">polyclonal</a> software package.
  </p>
</body>
</html>
EOF
        """

# Export the AA fitness heatmap dashboard as a WASM HTML file for GitHub Pages
rule export_dashboard:
    input:
        notebook="notebooks/aa_fitness_heatmap_dashboard.py",
        aa_fitness=f"{config['output_dir']}/aa_fitness_effects.csv",
        reference_aa=f"{config['output_dir']}/reference_aa.json"
    output:
        "docs/aa/index.html"
    log:
        f"{config['output_dir']}/logs/export_dashboard.log"
    shell:
        """
        marimo export html-wasm {input.notebook} -o docs/aa/ --mode run -f &> {log}
        sed -i 's/"theme": "light"/"theme": "dark"/g' docs/aa/index.html
        mkdir -p docs/aa/results
        cp {input.aa_fitness} docs/aa/results/
        cp {input.reference_aa} docs/aa/results/
        """

rule export_nt_dashboard:
    input:
        notebook="notebooks/nt_fitness_heatmap_dashboard.py",
        nt_fitness=f"{config['output_dir']}/nt_fitness_effects.csv",
        reference_nt=f"{config['output_dir']}/reference_nt.json",
    output:
        "docs/nt/index.html"
    log:
        f"{config['output_dir']}/logs/export_nt_dashboard.log"
    shell:
        """
        marimo export html-wasm {input.notebook} -o docs/nt/ --mode run -f &> {log}
        sed -i 's/"theme": "light"/"theme": "dark"/g' docs/nt/index.html
        mkdir -p docs/nt/results
        cp {input.nt_fitness} docs/nt/results/
        cp {input.reference_nt} docs/nt/results/
        """

# Compute mutation rates for subset trees (host, geographic, split-half)
rule compute_subset_rates:
    input:
        notebook="notebooks/compute_subset_rates.ipynb",
        mutation_counts=get_subset_mutation_count_files()
    output:
        subset_counts="{output_dir}/subset_counts.csv"
    log:
        "{output_dir}/logs/compute_subset_rates.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            compute_subset_rates.ipynb &> ../{log}
        """

# Compute fitness effects for subset trees using expected rates from global neutral model
rule compute_subset_fitness_effects:
    input:
        notebook="notebooks/compute_subset_fitness_effects.ipynb",
        subset_counts="{output_dir}/subset_counts.csv",
        expected_rates="{output_dir}/expected_rates.csv"
    output:
        subset_aa_fitness="{output_dir}/subset_aa_fitness_effects.csv"
    log:
        "{output_dir}/logs/compute_subset_fitness_effects.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            compute_subset_fitness_effects.ipynb &> ../{log}
        """

# Analyze fitness effects across subsets with scatter plots
rule analyze_subset_fitness_effects:
    input:
        notebook="notebooks/analyze_subset_fitness_effects.ipynb",
        subset_aa_fitness="{output_dir}/subset_aa_fitness_effects.csv"
    output:
        touch("{output_dir}/.analyze_subset_fitness_effects.done")
    log:
        "{output_dir}/logs/analyze_subset_fitness_effects.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            analyze_subset_fitness_effects.ipynb &> ../{log}
        """

# Render single-panel fitness-effect heatmap PNGs by re-using the same
# Altair specs that drive the marimo dashboards (aa/nt_fitness_heatmap_dashboard.py).
# Outputs go into results/figures/ alongside the other panel PNGs and feed
# compose_figures.ipynb.
HEATMAP_SNAPSHOT_NAMES = [
    "aa_heatmap",
    "nt_heatmap",
    "MP_5SS",
    "MP_3SS",
    "NS_5SS",
    "NS_3SS",
]


rule render_fitness_heatmaps:
    input:
        script="scripts/render_fitness_heatmaps.py",
        aa_csv="{output_dir}/aa_fitness_effects.csv",
        nt_csv="{output_dir}/nt_fitness_effects.csv",
        aa_ref="{output_dir}/reference_aa.json",
        nt_ref="{output_dir}/reference_nt.json",
    output:
        expand(
            "{{output_dir}}/figures/{name}.png",
            name=HEATMAP_SNAPSHOT_NAMES,
        )
    log:
        "{output_dir}/logs/render_fitness_heatmaps.log"
    shell:
        """
        python {input.script} \
            --aa-csv {input.aa_csv} --nt-csv {input.nt_csv} \
            --aa-ref {input.aa_ref} --nt-ref {input.nt_ref} \
            --output-dir {wildcards.output_dir}/figures \
            &> {log}
        """


# Compose multi-panel manuscript figures from individual panel PNGs.
# Reads panels produced by the upstream analyze_* notebooks, two manually
# staged PNGs in data/pngs/ (concept map, tree summary), the six fitness-effect
# heatmaps rendered by render_fitness_heatmaps, and the flu-usher leaves_per_tree.png.
rule compose_figures:
    input:
        notebook="notebooks/compose_figures.ipynb",
        genome_wide_done="{output_dir}/.analyze_genome_wide_rates.done",
        site_specific_done="{output_dir}/.analyze_site_specific_rates.done",
        fitness_effects_done="{output_dir}/.analyze_fitness_effects.done",
        leaves_per_tree=f"{config['data_dir']}/figures/leaves_per_tree.png",
        manual_pngs=expand(
            "data/pngs/{name}.png",
            name=["summary_concept_map", "tree_summary"],
        ),
        heatmap_snapshots=expand(
            "{{output_dir}}/figures/{name}.png",
            name=HEATMAP_SNAPSHOT_NAMES,
        ),
    output:
        touch("{output_dir}/.compose_figures.done")
    log:
        "{output_dir}/logs/compose_figures.log"
    shell:
        """
        mkdir -p {wildcards.output_dir}/figures/ms_figures && \
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            compose_figures.ipynb &> ../{log}
        """

# Summarize mutation filter statistics across all trees from log files
rule summarize_filter_logs:
    input:
        notebook="notebooks/summarize_filter_logs.ipynb",
        mutation_counts=get_all_mutation_count_files()
    output:
        touch("{output_dir}/.summarize_filter_logs.done")
    log:
        "{output_dir}/logs/summarize_filter_logs.log"
    shell:
        """
        cd notebooks && \
        jupyter nbconvert \
            --to notebook \
            --execute \
            --inplace \
            --ExecutePreprocessor.timeout=600 \
            summarize_filter_logs.ipynb &> ../{log}
        """