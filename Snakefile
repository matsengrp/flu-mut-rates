"""Top-level ``snakemake`` file that runs pipeline."""
import os
import textwrap
import yaml
from itertools import permutations
configfile: "config.yaml"



segs = config["files_to_grab"]["segments"].keys()
mut_types = list(map("".join,permutations("AGCT",2)))


rule all:
    input:       
        expand(
            #"results_{mat}/results/counts_by_mut_type_{seg}.png",
            "results_{mat}/results/tree_sizes_{seg}.csv",
            mat=config["files_to_grab"]["mat"],
            seg=segs
        )


rule get_flu_segments_from_web:
    """...get the files..."""
    output:
        expand(
            "results_{mat}/data/{seg_name}.{file}",
            mat=config["files_to_grab"]["mat"],
            seg_name=config["files_to_grab"]["segments"].keys(),
            file=["meta.tsv", "pb"]
        )
    params:
        the_urls=expand(
            "{base_url}{value}/{more}{value}{file}", 
            base_url=config["files_to_grab"]["base_url"],
            value=config["files_to_grab"]["segments"].values(),
            more=config["files_to_grab"]["more_url"],
            file=[
                config["files_to_grab"]["meta"],
                config["files_to_grab"]["pb"]
            ]
        )
    run:
        for url, dest in zip(params.the_urls, output):
            shell(f"wget -O - {url} | gunzip -c > {dest}")


metas = expand("{seg}.meta.tsv",seg=segs)
rule find_roots:
    """... Determine new root nodes. """
    input:
        ["results_{mat}/data/" + file for file in metas]
    output:
        path="results_{mat}/data/root_ids.csv"
    script:
        "scripts/find_roots.py"


pbs = expand("{seg}.pb",seg=segs)
rule reroot_trees:
    """...reroot to common old node..."""
    input:
        tree_paths=["results_{mat}/data/" + file for file in pbs],
        root_ids="results_{mat}/data/root_ids.csv"
    output:
        ["results_{mat}/data/rerooted_" + file for file in pbs]
    conda:
        "flu-syn-rates-usher"
    shell:
        """
        # Create a bash array from the input files
        LIST1=({input.tree_paths})    
        LIST2=({output})
        COUNTER=0

        while IFS=',' read -r c1 c2; do    
            inpath=${{LIST1[${{COUNTER}}]}}
            outpath=${{LIST2[${{COUNTER}}]}}                
            matUtils extract -i $inpath --reroot "$c2" -o $outpath
            #matUtils reroot -i $inpath -o $outpath -k "$c2"--preserve-branch-length       
            let COUNTER=1+$COUNTER
        done < "{input.root_ids}"
        """


coding_site_files = expand('coding_sites_{seg}.csv', seg=segs)
fasta_files = expand('inferred_{seg}.fasta', seg=segs)
gtf_files = expand('ref_{seg}.gtf', seg=segs)
rule makes_coding_dicts:
    input:
        original_tree_paths=["results_{mat}/data/" + file for file in pbs],
        rerooted_tree_paths=["results_{mat}/data/rerooted_" + file for file in pbs],
        root_ids="results_{mat}/data/root_ids.csv"
    output:
        coding_site_paths=["results_{mat}/data/"+file for file in coding_site_files],
        fasta_paths=["results_{mat}/data/"+file for file in fasta_files],
        gtf_paths=["results_{mat}/data/"+file for file in gtf_files]
    script:
        "scripts/make_coding_sites.py"


bad_site_files = expand('bad_sites_{seg}.csv', seg=segs)
secondary_structs_files = expand('secondary_struct_{seg}.csv', seg=segs)
rates_tables_files = expand('rates_{seg}.csv', seg=segs)
rule make_placeholder_files:
    output: 
        bad_sites_paths=["results_{mat}/data/"+file for file in bad_site_files],
        secondary_structs=["results_{mat}/data/"+file for file in secondary_structs_files], 
        rates_tables=["results_{mat}/data/"+file for file in rates_tables_files],
    script:
        "scripts/make_placeholder_files.py"


all_counts_files = expand("all_counts_{seg}.csv", seg=segs)
all_pcps_files = expand("all_pcps_{seg}.csv", seg=segs)
nuc_conservation_files = expand("nt_conservation_{seg}.csv", seg=segs)
codon_conservation_files = expand("codon_conservation_{seg}.csv", seg=segs)
tree_size_files = expand("tree_sizes_{seg}.csv", seg=segs)
rule make_count_dfs:
    input:
        rerooted_tree_paths=["results_{mat}/data/rerooted_" + file for file in pbs],
        coding_site_paths=["results_{mat}/data/"+file for file in coding_site_files],
        fasta_paths=["results_{mat}/data/"+file for file in fasta_files],
        gtf_paths=["results_{mat}/data/"+file for file in gtf_files],
        bad_sites_paths=["results_{mat}/data/"+file for file in bad_site_files],
        secondary_structs=["results_{mat}/data/"+file for file in secondary_structs_files], 
        rates_tables=["results_{mat}/data/"+file for file in rates_tables_files],
    output:
        all_counts_paths=["results_{mat}/results/"+file for file in all_counts_files],
        all_pcps_paths=["results_{mat}/results/"+file for file in all_pcps_files],
        nuc_conservation_paths=["results_{mat}/results/"+file for file in nuc_conservation_files],
        codon_conservation_paths=["results_{mat}/results/"+file for file in codon_conservation_files],
        tree_size_paths=["results_{mat}/results/"+file for file in tree_size_files]
    script:
        "scripts/make_count_dfs.py"


boxplot_files = expand("counts_by_mut_type_{seg}.png", seg=segs)
mut_type_dirs = expand("scatter_plots_{seg}/", seg=segs)
mut_type_files = expand("scatter_plots_{seg}/{mt}.png", seg=segs, mt=mut_types)
rule make_basic_plots:
    input:
        all_counts_paths=["results_{mat}/results/"+file for file in all_counts_files]
    output:
        boxplot_paths = ["results_{mat}/results/"+file for file in boxplot_files],
        mut_type_paths = ["results_{mat}/results/"+file for file in mut_type_files]
    params:
        mut_type_dirs = ["results_{mat}/results/"+file for file in mut_type_dirs]
    script:
        "scripts/basic_plots.py"