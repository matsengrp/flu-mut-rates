import pandas as pd
import math
from collections import defaultdict
import bte
from ExpectedCalc import apply_muts
import requests

# from Bio import Entrez, SeqIO


old_trees = snakemake.input.original_tree_paths
new_trees = snakemake.input.rerooted_tree_paths
coding_site_paths = snakemake.output.coding_site_paths
fasta_paths = snakemake.output.fasta_paths
gtf_paths = snakemake.output.gtf_paths
roots = snakemake.input.root_ids


def get_sequence(seq_id, email="AnOtherExample%40example.com"):
    ### Preferred approach using their API.
    ### An email address is necessary, but it need not be valid.
    # Entrez.email = "A.N.Other@example.com"
    # wrapper = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    # the_seq = str(SeqIO.read(wrapper, "fasta").seq)
    # wrapper.close()

    ###But that has a certificate problem, so we do it manually.
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    url += f"{seq_id}&rettype=fasta&retmode=text&tool=biopython&email={email}"
    response = requests.get(url)
    the_seq = "".join(response.text.split("\n")[1:])
    return the_seq


def get_valid_codon_pos(seq_id, email="AnOtherExample%40example.com"):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    url += f"{seq_id}&rettype=ft&retmode=text&tool=biopython&email={email}"
    response = requests.get(url)
    start, stop = next(
        map(int, line.split("\t")[:2])
        for line in map(str.strip, response.text.split("\n")[1:])
        if line.endswith("gene")
    )
    return start, stop


def guess_old_root_seq_id(tree):
    seq_id = next(
        x[-2]
        for c in tree.root.children
        if len(c.mutations) == 0 and len((x := c.id.split("|"))) >= 3
    )
    return seq_id


def gather_and_write(
    old_tree, new_tree, gene, root_id, coding_site_out, fasta_out, gtf_out
):
    print(f"Processing segment {gene}")

    o_tree = bte.MATree(old_tree)
    old_root_seq_id = guess_old_root_seq_id(o_tree)
    o_ref_seq = get_sequence(old_root_seq_id)
    gene_start, gene_end = get_valid_codon_pos(old_root_seq_id)
    with open("old_ref_seq.fasta", "w") as the_file:
        the_file.write(o_ref_seq)

    r_tree = bte.MATree(new_tree)
    new_root_seq_id = root_id.split("|")[-2]
    r_ref_seq = get_sequence(new_root_seq_id)
    muts = {int(mut[1:-1]): mut[-1] for mut in set(o_tree.get_haplotype(root_id))}

    ref_seq = o_ref_seq if len(muts) == 0 else apply_muts(o_ref_seq, muts)
    # The inferred ref_seq should match the r_ref_seq once aligned and whatever else.

    coding_sites_dict = defaultdict(list)
    for site in range(1, len(ref_seq) + 1):
        coding_sites_dict["site"].append(site)
        if gene_start <= site <= gene_end:
            offset = site - gene_start + 1
            rem = offset % 3
            coding_sites_dict["codon_position"].append(3 if rem == 0 else rem)
            coding_sites_dict["codon_site"].append(math.ceil(offset / 3))
            coding_sites_dict["gene"].append(gene)
        else:
            coding_sites_dict["codon_position"].append("noncoding")
            coding_sites_dict["codon_site"].append("noncoding")
            coding_sites_dict["gene"].append("noncoding")

    coding_sites_df = pd.DataFrame(coding_sites_dict)
    coding_sites_df.to_csv(coding_site_out)
    with open(fasta_out, "w") as the_file:
        the_file.write(">\n")
        the_file.write(ref_seq)
    gtf_string = f"{new_root_seq_id}\tncbiGenes.genePred\tCDS\t{gene_start}\t{gene_end}"
    gtf_string += f'\t.\t+\t0\tgene_id "{gene}";'
    with open(gtf_out, "w") as the_file:
        the_file.write(gtf_string)

    return None


# Call those functions.
with open(roots) as the_file:
    root_dict = {
        (x := line.split(","))[0]: x[1] for line in map(str.strip, the_file.readlines())
    }
params_list = zip(
    old_trees,
    new_trees,
    root_dict.keys(),
    root_dict.values(),
    coding_site_paths,
    fasta_paths,
    gtf_paths,
)
any((gather_and_write(*p) for p in params_list))
