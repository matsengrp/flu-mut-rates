# Load configuration
configfile: "config/config.yaml"

# Helper function to get subtypes for each segment
def get_subtypes_for_segment(segment):
    if segment == "HA":
        return config["ha_subtypes"]
    elif segment == "NA":
        return config["na_subtypes"]
    else:
        return ["all"]

# Generate all segment-subtype-host combinations
segment_subtype_host_combinations = []
for segment in config["segments"]:
    for subtype in get_subtypes_for_segment(segment):
        for host in config["host_groups"]:
            segment_subtype_host_combinations.append((segment, subtype, host))

# Final output files
final_outputs = []
for segment, subtype, host in segment_subtype_host_combinations:
    final_outputs.extend([
        f"{config['output_dir']}/{segment}/{subtype}/{host}/mutation_counts.csv",
        f"{config['output_dir']}/{segment}/{subtype}/{host}/parent_child_pairs.csv"
    ])

# Add aligned proteins outputs (only for HA and NA segments)
final_outputs.extend(expand(
    "{output_dir}/aligned_proteins/{segment}",
    output_dir=config["output_dir"],
    segment=["HA", "NA"]
))

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
        "logs/{segment}_{subtype}_coding_sites.log"
    shell:
        """
        python scripts/make_coding_sites.py \
            --ref_fasta {input.ref_fasta} \
            --gff_file {input.gff_file} \
            --subtype {wildcards.subtype} \
            --segment {wildcards.segment} \
            --output_directory {wildcards.output_dir}/{wildcards.segment}/{wildcards.subtype} &> {log}
        """

# Count mutations along tree
rule count_mutations:
    input:
        tree_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/host_specific_trees/{wildcards.host}_tree.pb.gz",
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_root.fasta",
        gtf_path=lambda wildcards: f"{config['data_dir']}/{wildcards.segment}/{wildcards.subtype}/curated_reference.gtf"
    output:
        all_counts_path="{output_dir}/{segment}/{subtype}/{host}/mutation_counts.csv",
        all_pcps_path="{output_dir}/{segment}/{subtype}/{host}/parent_child_pairs.csv"
    log:
        "logs/{segment}_{subtype}_{host}_mutation_counts.log"
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
        "logs/align_proteins_{segment}.log"
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