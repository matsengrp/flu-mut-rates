# import pandas as pd
from pandas import DataFrame as DF
from collections import defaultdict, Counter
import bte
from ExpectedCalc import PossibleMutations, apply_muts

# from MATWrapper import apply_muts


rerooted_tree_paths = snakemake.input.rerooted_tree_paths
coding_site_paths = snakemake.input.coding_site_paths
fasta_paths = snakemake.input.fasta_paths
gtf_paths = snakemake.input.gtf_paths
bad_sites_paths = snakemake.input.bad_sites_paths
secondary_structs = snakemake.input.secondary_structs
rates_tables = snakemake.input.rates_tables
nuc_count_paths = snakemake.output.nuc_count_paths
nuc_conservation_paths = snakemake.output.nuc_conservation_paths
codon_count_paths = snakemake.output.codon_count_paths
codon_conservation_paths = snakemake.output.codon_conservation_paths
possible_muts_paths = snakemake.output.possible_muts_paths
all_count_paths = snakemake.output.all_count_paths


def get_3mer_seq_context(nt_site, seq):
    return seq[nt_site - 2 : nt_site + 1]


def fails_filters(node, ref_seq, nt_muts, codon_muts):
    fail = node.parent is None
    fail |= len(nt_muts) > 4
    fail |= sum(ref_seq[int(mut[1:-1]) - 1] == mut[-1] for mut in nt_muts) > 1
    fail |= len({(mut.gene, mut.aa_index) for mut in codon_muts}) < len(codon_muts)
    return fail


def write_counts(
    tree_path,
    coding_site_path,
    fasta_path,
    gtf_path,
    bad_sites,
    secondary,
    rates,
    nuc_count_path,
    nuc_conserv_path,
    codon_count_path,
    codon_conserv_path,
    possible_muts_path,
    all_count_path,
):
    possible_muts = PossibleMutations(
        fasta_path, coding_site_path, bad_sites, secondary, rates, False
    )
    ref_seq = possible_muts.ref_seq
    ref_df = possible_muts.possible_mutations_df()
    ref_df.rename(columns={"codon_sites": "codon_site"}, inplace=True)
    ref_df.drop(columns=["ID"], inplace=True)
    ref_df["gene_codon_site_id"] = ref_df["gene"] + "-" + ref_df["codon_site"]
    ref_df["seq_context"] = ref_df.site.apply(get_3mer_seq_context, args=(ref_seq,))
    ref_df.to_csv(possible_muts_path)

    tree = bte.MATree(tree_path)
    translations = tree.translate(gtf_path, fasta_path)

    # Cycle over nodes and record codon mutation counts
    nodes = tree.depth_first_expansion()
    print("Number of nodes in tree:", len(nodes))
    nt_mut_counts = defaultdict(int)
    codon_mut_counts = defaultdict(int)
    conserv_nt_ctr = Counter()
    conserv_codon_ctr = Counter()
    n_nodes_passing_filters = 0
    n_internal_nodes_analyzed = 0
    for i, node in enumerate(nodes):
        # Get lists of nt-level and codon-level mutations
        nt_muts = node.mutations
        codon_muts = translations[node.id] if node.id in translations else []
        if fails_filters(node, ref_seq, nt_muts, codon_muts):
            continue

        n_nodes_passing_filters += 1

        # Reconstruct the parent node's genome for use in getting each mutation's
        # sequence context
        node_haplotype = tree.get_haplotype(node.id)
        parent_node_haplotype = tree.get_haplotype(node.parent.id)
        parent_muts = {int(mut[1:-1]): mut[-1] for mut in parent_node_haplotype}
        parent_node_seq = apply_muts(ref_seq, parent_muts)

        # Record counts of nt-level mutations
        for nt_mut in nt_muts:
            seq_context = get_3mer_seq_context(int(nt_mut[1:-1]), parent_node_seq)
            ref_seq_context = get_3mer_seq_context(int(nt_mut[1:-1]), ref_seq)
            nt_mut_id = f"{nt_mut}-{seq_context}-{ref_seq_context}"
            nt_mut_counts[nt_mut_id] += 1

        # Record counts of codon-level mutations
        for codon_mut in codon_muts:
            nt_mut = codon_mut.nuc
            seq_context = get_3mer_seq_context(int(nt_mut[1:-1]), parent_node_seq)
            ref_seq_context = get_3mer_seq_context(int(nt_mut[1:-1]), ref_seq)
            codon_mut_id = f"{codon_mut.gene}-{codon_mut.aa_index}-{codon_mut.original_codon}-{codon_mut.alternative_codon}-{seq_context}-{ref_seq_context}"
            codon_mut_counts[codon_mut_id] += 1

        # Record rolling counts of all nt-level and codon-level mutations in a node
        # relative to the reference sequence in order to determine the level of
        # conservation at each site. Only do this for internal nodes.
        if not node.is_leaf():
            n_internal_nodes_analyzed += 1
            node_counts = Counter([mut[1:-1] for mut in node_haplotype])
            conserv_nt_ctr += node_counts
            q = "nuc.isin(@node_haplotype)"
            node_counts = Counter(set(ref_df.query(q).gene_codon_site_id))
            conserv_codon_ctr += node_counts

        if i % 10000 == 0:
            print(i)

    print("Number of nodes that passed filters:", n_nodes_passing_filters)
    print("Number of internal nodes analyzed:", n_internal_nodes_analyzed)

    # Store mutation counts in dataframes and write out.
    cols = ["nt_mut_id", "counts"]
    nt_mut_counts_df = DF.from_records(list(nt_mut_counts.items()), columns=cols)
    new_cols = ["nt_mut", "mut_seq_context", "ref_seq_context"]
    new_vals = nt_mut_counts_df["nt_mut_id"].str.extract("(\w+)-(\w+)-(\w+)")
    nt_mut_counts_df[new_cols] = new_vals
    nt_mut_counts_df.drop(columns=["nt_mut_id"], inplace=True)
    nt_mut_counts_df.to_csv(nuc_count_path, index=False)

    cols = ["codon_mut_id", "counts"]
    codon_mut_counts_df = DF.from_records(list(codon_mut_counts.items()), columns=cols)
    new_cols = ["gene", "codon_site", "wt_codon", "mut_codon", "mut_seq_context"]
    new_cols.append("ref_seq_context")
    id_format = "(\w+)-(\d+)-(\w+)-(\w+)-(\w+)-(\w+)"
    new_vals = codon_mut_counts_df["codon_mut_id"].str.extract(id_format)
    codon_mut_counts_df[new_cols] = new_vals
    codon_mut_counts_df.drop(columns=["codon_mut_id"], inplace=True)
    codon_mut_counts_df.to_csv(codon_count_path, index=False)

    join_cols = ["codon_site", "wt_codon", "mut_codon", "gene"]
    all_counts_df = ref_df.merge(codon_mut_counts_df, on=join_cols, how="left")
    all_counts_df["mut_type"] = all_counts_df.wt_nt + all_counts_df.mut_nt
    all_counts_df["codon_site"] = all_counts_df.codon_site.astype(int)
    all_counts_df["counts"] = all_counts_df.counts.fillna(0).astype(int)
    all_counts_df.to_csv(all_count_path, index=False)

    pairs = (conserv_nt_ctr, nuc_conserv_path), (conserv_codon_ctr, codon_conserv_path)
    for ctr, path in pairs:
        the_df = DF.from_dict(ctr, orient="index")
        the_df.reset_index(inplace=True)
        the_df.rename(columns={"index": "site", 0: "count"}, inplace=True)
        the_df["frac_mut"] = the_df["count"] / n_internal_nodes_analyzed
        the_df["frac_cons"] = 1 - the_df["frac_mut"]
        the_df.sort_values("count", ascending=False, inplace=True)
        the_df.to_csv(path, index=False)

    return None


params = zip(
    rerooted_tree_paths,
    coding_site_paths,
    fasta_paths,
    gtf_paths,
    bad_sites_paths,
    secondary_structs,
    rates_tables,
    nuc_count_paths,
    nuc_conservation_paths,
    codon_count_paths,
    codon_conservation_paths,
    possible_muts_paths,
    all_count_paths,
)
any(write_counts(*p) for p in params)
