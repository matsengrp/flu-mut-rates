# Load configuration
configfile: "config/config.yaml"

# Final output files
final_outputs = expand(
    ["{output_dir}/{subtype}/{segment}/mutation_counts.csv",
     "{output_dir}/{subtype}/{segment}/parent_child_pairs.csv"],
    output_dir=config["output_dir"],
    subtype=config["subtypes"],
    segment=config["segments"]
)

# Main rule to define target outputs
rule all:
    input:
        final_outputs

# Create coding sites file
rule make_coding_sites:
    input:
        ref_fasta=lambda wildcards: config["reference_format"].format(
            data_dir=config["data_dir"],
            subtype=wildcards.subtype,
            segment=wildcards.segment
        ),
        gff_file=lambda wildcards: config["gff_format"].format(
            data_dir=config["data_dir"],
            subtype=wildcards.subtype,
            segment=wildcards.segment
        )
    output:
        coding_sites="{output_dir}/{subtype}/{segment}/coding_sites.csv"
    log:
        "{output_dir}/logs/{subtype}_{segment}_coding_sites.log"
    shell:
        """
        python scripts/make_coding_sites.py \
            --ref_fasta {input.ref_fasta} \
            --gff_file {input.gff_file} \
            --output {output.coding_sites} 2> {log}
        """

# Count mutations along tree
rule count_mutations:
    input:
        tree_path=lambda wildcards: config["tree_format"].format(
            data_dir=config["data_dir"],
            subtype=wildcards.subtype,
            segment=wildcards.segment
        ),
        coding_site_path=rules.make_coding_sites.output.coding_sites,
        fasta_path=lambda wildcards: config["reference_format"].format(
            data_dir=config["data_dir"],
            subtype=wildcards.subtype,
            segment=wildcards.segment
        ),
        gtf_path=lambda wildcards: config["gtf_format"].format(
            data_dir=config["data_dir"],
            subtype=wildcards.subtype,
            segment=wildcards.segment
        )
    output:
        all_counts_path="{output_dir}/{subtype}/{segment}/mutation_counts.csv",
        all_pcps_path="{output_dir}/{subtype}/{segment}/parent_child_pairs.csv"
    log:
        "{output_dir}/logs/{subtype}_{segment}_mutation_counts.log"
    shell:
        """
        python scripts/make_count_dfs.py \
            --tree_path {input.tree_path} \
            --coding_site_path {input.coding_site_path} \
            --fasta_path {input.fasta_path} \
            --gtf_path {input.gtf_path} \
            --all_counts_path {output.all_counts_path} \
            --all_pcps_path {output.all_pcps_path} 2> {log}
        """