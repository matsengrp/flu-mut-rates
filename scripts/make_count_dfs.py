from pandas import DataFrame as DF
from collections import defaultdict, Counter
import bte
from ExpectedCalc import PossibleMutations, apply_muts
from time import time
from multiprocessing import Pool
from math import ceil

# from MATWrapper import apply_muts


class CountsHelper:
    """
    Helper class for calculating and writing per site mutation counts and conservation.

    Attributes:
        int_nodes (list): The internal nodes of the tree.
        key_names (tuple): Column names that also serve as keys for Counters.
        n_int_nodes (int): The number of internal nodes of the tree.
        n_nodes (int): The number of nodes in the tree.
        nodes (list): The nodes of the tree.
        poss_muts (PossibleMutations): Possible Mutations instance for the data.
        ref_df (pd.DataFrame) = DataFrame of possible mutations from the reference
            sequence, along with additonal information.
        ref_seq (str): The reference sequence.
        translations (defaultdict): A default dictionary mapping a node id to the codon
            mutations at the node.
        tree (bte.MATree): Instance of the tree.
    """

    def __init__(
        self,
        tree_path,
        coding_site_path,
        fasta_path,
        gtf_path,
        bad_sites,
        secondary,
        rates,
    ):
        """
        Parameters:
            tree_path (str): The tree.
            coding_site_path (str): The coding sites dataframe (nucleotide site to codon
                and gene).
            fasta_path (str): The fasta for the reference sequence.
            gtf_path (str):
            bad_sites (str): The file listing bad sites for initializing self.poss_muts.
                Currently this is a dummy file, no known bad sites.
            secondary (str): The file listing secondary structure predictions for
                initializing self.poss_muts. Currently this is a dummy file, no known
                secondary structure.
            rates (str): The file of rates for different mutation types for initializing
                self.poss_muts. Currently this is a dummy file, no master rates table.
        """

        self.poss_muts = PossibleMutations(
            fasta_path, coding_site_path, bad_sites, secondary, rates, False
        )
        self.ref_seq = self.poss_muts.ref_seq

        ref_df = self.poss_muts.possible_mutations_df()
        self.key_names = (*ref_df.columns, "parent_motif")
        ref_df["gene_codon_site_id"] = ref_df["gene"] + "-" + ref_df["codon_sites"]
        ref_df["seq_context"] = ref_df.site.apply(self.ref_motif)
        self.ref_df = ref_df

        self.tree = bte.MATree(tree_path)
        self.translations = defaultdict(list, self.tree.translate(gtf_path, fasta_path))
        self.nodes = self.tree.depth_first_expansion()
        self.int_nodes = [n for n in self.nodes if not n.is_leaf()]
        self.n_nodes = len(self.nodes)
        self.n_int_nodes = len(self.int_nodes)

    def ref_motif(self, s):
        """Return the 3-mer motif centered at the site. Site indices begin at one."""
        return self.ref_seq[s - 2 : s + 1]

    def fails_filters(self, nt_muts, codon_muts):
        """
        Returns true if the mutations at a node fail the filter. The filter requires
        that there are between one and four mutations; there is at most one reversion to
        the reference sequence; and no two mutations target the same codon.

        Parameters:
            nt_mut (list): List of nucleotide mutations, string entries like "A1234G".
            codon_muts (list): List of codon mutations, entries are bte.AAChange
                objects.
        """
        fail = len(nt_muts) > 4
        fail |= len(nt_muts) == 0
        fail |= sum(self.ref_seq[int(mut[1:-1]) - 1] == mut[-1] for mut in nt_muts) > 1
        fail |= len({(mut.gene, mut.aa_index) for mut in codon_muts}) < len(codon_muts)
        return fail

    def count_mutations_from_parent(self, parent):
        """
        Count mutations from the parent node to all of its child nodes. Return four
        values, which are the number of child nodes that pass the filter;
        a Counter mapping a tuple of values for self.key_names to the number of times
        the associating mutation type occured; a Counter mapping a tuple of values for
        self.key_names to the sum of branch lengths where the mutation type is possible;
        and a list of parent child pairs, where each entry is a triple with of the
        parent node id, the parent node sequence, and a list of pairs of child node id
        and child node sequence.
        """
        parent_node_haplotype = self.tree.get_haplotype(parent.id)
        p_muts = {int(mut[1:-1]): mut[-1] for mut in parent_node_haplotype}

        # If speed is an issue, we can rewrite the possible_muts_from_founder_muts
        # method to avoid dataframes.
        counts_df, parent_seq = self.poss_muts.possible_muts_from_founder_muts(p_muts)
        parent_motif = lambda s: parent_seq[s - 2 : s + 1]
        counts_df["parent_motif"] = counts_df.site.apply(parent_motif)
        counts_df["actual_count"] = 0
        counts_df["branch_length"] = 0
        pcps = (parent.id, parent_seq, [])

        n_filtered = 0
        for node in parent.children:
            nt_muts = node.mutations
            codon_muts = self.translations[node.id]
            if self.fails_filters(nt_muts, codon_muts):
                continue

            n_filtered += 1
            pcps[2].append((node.id, nt_muts))
            counts_df["branch_length"] += len(nt_muts)
            to_increment = counts_df.query("nuc.isin(@node.mutations)").index
            counts_df.loc[to_increment, "actual_count"] += 1

        if n_filtered != 0:
            records = counts_df.to_records(index=False).tolist()
            counts_dict = Counter({x[:-2]: x[-2] for x in records})
            branches_dict = Counter({x[:-2]: x[-1] for x in records})
        else:
            counts_dict, branches_dict = Counter(), Counter()

        return n_filtered, counts_dict, branches_dict, pcps

    def count_sites_on_node(self, node):
        """
        Count the sites targetted by mutations on this node. Return two Counters, the
        first maps a nucleotide site to the number of mutations at that site at this
        node; the second maps a string identifier of the form "{gene}-{codon_site}" to
        the number of mutations at that codon site at this node.
        """
        node_haplotype = self.tree.get_haplotype(node.id)
        nt_site_count = Counter([mut[1:-1] for mut in node_haplotype])
        view = self.ref_df.query("nuc.isin(@node_haplotype)")
        codon_site_count = Counter(set(view.gene_codon_site_id))
        return nt_site_count, codon_site_count

    def mut_counters_to_df(self, actual_counts, branch_lengths):
        """
        Return a DataFrame with data from two counters like those returned by
        self.count_mutations_from_parent.
        """
        zipped = zip(actual_counts.items(), branch_lengths.values())
        data = [(*k, v, w) for (k, v), w in zipped]
        columns = [*self.key_names, "actual_count", "branch_length"]
        the_df = DF(data, columns=columns)
        return the_df

    def site_counter_to_df(self, counter, n_filtered):
        """
        Return a dataframe with site conservation information from a Counter.
        """
        the_df = DF.from_dict(counter, orient="index")
        the_df.reset_index(inplace=True)
        the_df.rename(columns={"index": "site", 0: "mut_count"}, inplace=True)
        the_df["frac_mut"] = the_df.mut_count / n_filtered
        the_df["frac_cons"] = 1 - the_df.frac_mut
        the_df.sort_values("mut_count", ascending=False, inplace=True)
        return the_df

    def compute_conservation_on_tree(self):
        """
        Compute the site conservation along the tree. Return the number of internal
        nodes passing the filter, a dataframe with nucleotide site conversation, and a
        dataframe with codon site conservation.
        """
        # Use counters rather than a dataframe, since summing together many.
        n_filtered, nt_site_counter, codon_site_counter = 0, Counter(), Counter()
        t0 = -time()
        for i, node in enumerate(self.int_nodes, 1):
            if i % 1000 == 0:
                print(f"processing node {i}; {t0+time():0.1f} seconds")

            nt_muts = node.mutations
            codon_muts = self.translations[node.id]
            if self.fails_filters(nt_muts, codon_muts):
                continue
            n_filtered += 1

            nt_site_count, codon_site_count = self.count_sites_on_node(node)
            nt_site_counter.update(nt_site_count)
            codon_site_counter.update(codon_site_count)

        nt_site_df = self.site_counter_to_df(nt_site_counter, n_filtered)
        codon_site_df = self.site_counter_to_df(codon_site_counter, n_filtered)
        return n_filtered, nt_site_df, codon_site_df

    def count_mutations_on_tree(self):
        """
        Count the number of mutations along the tree. Return the number of nodes passing
        the filter, a dataframe of per site nucleotide mutations, and a dataframe of
        parent child pairs.
        """
        # Use counters rather than a dataframe, since summing together many.
        n_filtered, all_actual_counts, all_branch_lengths = 0, Counter(), Counter()
        all_pcps = []
        t0 = -time()
        for i, node in enumerate(self.int_nodes, 1):
            if i % 1000 == 0:
                print(f"processing node {i}; {t0+time():0.1f} seconds")

            count, actual_counts, branch_lens, pcps = self.count_mutations_from_parent(
                node
            )
            n_filtered += count
            all_actual_counts.update(actual_counts)
            all_branch_lengths.update(branch_lens)
            if count != 0:
                all_pcps.append(pcps)

        all_counts_df = self.mut_counters_to_df(all_actual_counts, all_branch_lengths)
        all_pcps_df = self.pcp_list_to_df(all_pcps)
        return n_filtered, all_counts_df, all_pcps_df

    def parallel_count_mutations_on_tree(self, processes=8, batch_size=1000):
        """
        ...would be nice, but not there yet...
        """
        raise NotImplementedError(
            "Can't multiprocess with bte tree, need to do extra book keeping."
        )
        n_filtered, all_actual_counts, all_branch_lengths = 0, Counter(), Counter()
        batch_count = ceil(self.n_int_nodes / batch_size)
        t0 = -time()
        with Pool(processes=processes) as pool:
            for b in range(batch_count):
                if b == batch_count - 1:
                    nodes = self.int_nodes[b * batch_size : (b + 1) * batch_size]
                else:
                    nodes = self.int_nodes[b * batch_size :]
                batch_result = pool.map(self.count_mutations_from_parent, nodes)
                for count, actual_counts, branch_lens in batch_result:
                    n_filtered += count
                    all_actual_counts.update(actual_counts)
                    all_branch_lengths.update(branch_lens)

            print(
                f"{b} batches of {batch_size} nodes processed in {t0+time():0.1f} seconds"
            )
        all_counts_df = self.mut_counters_to_df(all_actual_counts, all_branch_lengths)
        return n_filtered, all_counts_df

    def pcp_list_to_df(self, pcps):
        """
        Turn a list of pcps, like those returned by self.count_mutations_from_parent,
        into a DataFrame.
        """
        mut_dict = lambda muts: {int(mut[1:-1]): mut[-1] for mut in muts}
        child_seq = lambda seq, muts: apply_muts(seq, mut_dict(muts))
        cols = ["parent_name", "child_name", "parent", "child", "branch_length"]
        rows = (
            (p_id, c_id, p_seq, child_seq(p_seq, c_muts), len(c_muts))
            for (p_id, p_seq, child_entries) in pcps
            for c_id, c_muts in child_entries
        )
        the_df = DF(rows, columns=cols)
        return the_df

    def make_and_write_all_counts(
        self,
        all_counts_path,
        all_pcps_path,
        nuc_conserv_path,
        codon_conserv_path,
        tree_size_path,
    ):
        """Calculate various counts and write to file."""

        print(f"Number of nodes in tree: {self.n_nodes}")
        print(f"Number of internal nodes in tree: {self.n_int_nodes}")

        print(f"Counting mutations on tree...")
        n_filtered, all_counts_df, all_pcps_df = self.count_mutations_on_tree()
        print(f"Number of nodes passing filter: {n_filtered}")

        print(f"Calculating site conservation on tree...")
        n_int_filtered, nt_site_df, codon_site_df = self.compute_conservation_on_tree()
        print(f"Number of internal nodes passing filter: {n_filtered}")

        all_counts_df.to_csv(all_counts_path)
        all_pcps_df.to_csv(all_pcps_path)
        nt_site_df.to_csv(nuc_conserv_path)
        codon_site_df.to_csv(codon_conserv_path)

        with open(tree_size_path, "w") as the_file:
            header = "nodes,internal_nodes,filtered_nodes,filtered_internal_nodes\n"
            row = f"{self.n_nodes},{self.n_int_nodes},{n_filtered},{n_int_filtered}\n"
            the_file.write(header)
            the_file.write(row)

        return None


def write_counts(
    tree_path,
    coding_site_path,
    fasta_path,
    gtf_path,
    bad_sites,
    secondary,
    rates,
    all_counts_path,
    all_pcps_path,
    nuc_conserv_path,
    codon_conserv_path,
    tree_size_path,
):
    p1 = tree_path, coding_site_path, fasta_path, gtf_path, bad_sites, secondary, rates
    p2 = (
        all_counts_path,
        all_pcps_path,
        nuc_conserv_path,
        codon_conserv_path,
        tree_size_path,
    )
    CountsHelper(*p1).make_and_write_all_counts(*p2)
    return None


if __name__ == "__main__":
    # In the main hook, so the rest of the file can be imported by other files.

    rerooted_tree_paths = snakemake.input.rerooted_tree_paths
    coding_site_paths = snakemake.input.coding_site_paths
    fasta_paths = snakemake.input.fasta_paths
    gtf_paths = snakemake.input.gtf_paths
    bad_sites_paths = snakemake.input.bad_sites_paths
    secondary_structs = snakemake.input.secondary_structs
    rates_tables = snakemake.input.rates_tables
    all_counts_paths = snakemake.output.all_counts_paths
    all_pcps_paths = snakemake.output.all_pcps_paths
    nuc_conservation_paths = snakemake.output.nuc_conservation_paths
    codon_conservation_paths = snakemake.output.codon_conservation_paths
    tree_size_paths = snakemake.output.tree_size_paths

    params = zip(
        rerooted_tree_paths,
        coding_site_paths,
        fasta_paths,
        gtf_paths,
        bad_sites_paths,
        secondary_structs,
        rates_tables,
        all_counts_paths,
        all_pcps_paths,
        nuc_conservation_paths,
        codon_conservation_paths,
        tree_size_paths,
    )
    any(write_counts(*p) for p in params)
