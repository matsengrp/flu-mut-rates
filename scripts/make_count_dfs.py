import pandas as pd
from collections import defaultdict, Counter
import bte
from ExpectedCalc import PossibleMutations, apply_muts
from time import time
from multiprocessing import Pool
from math import ceil
import argparse

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
        gtf_path
    ):
        """
        Parameters:
            tree_path (str): The tree.
            coding_site_path (str): The coding sites dataframe (nucleotide site to codon
                and gene).
            fasta_path (str): The fasta for the reference sequence.
            gtf_path (str):
        """

        # Given the reference sequence and the a table of coding sites, initialize a
        # PossibleMutations instance for computing a dataframe of possible mutations.
        self.poss_muts = PossibleMutations(fasta_path, coding_site_path)
        self.ref_seq = self.poss_muts.ref_seq

        # Make a dataframe of all possible mutations to the reference sequence
        ref_df = self.poss_muts.possible_mutations_df()
        self.key_names = (*ref_df.columns, "parent_motif")
        ref_df["gene_codon_site_id"] = ref_df["gene"] + "-" + ref_df["codon_site"].astype(str)
        self.ref_df = ref_df

        # Read in the tree, translate mutations, and make a list of nodes
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
        that there are between one and four mutations and no two mutations target the same codon.

        Parameters:
            nt_mut (list): List of nucleotide mutations, string entries like "A1234G".
            codon_muts (list): List of codon mutations, entries are bte.AAChange
                objects.
        """
        fail = len(nt_muts) > 4
        fail |= len(nt_muts) == 0
        # fail |= sum(self.ref_seq[int(mut[1:-1]) - 1] == mut[-1] for mut in nt_muts) > 1
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

        # Get the set of mutations in the partent node relative to the reference
        parent_node_haplotype = self.tree.get_haplotype(parent.id)
        p_muts = {int(mut[1:-1]): mut[-1] for mut in parent_node_haplotype}

        # If speed is an issue, we can rewrite the possible_muts_from_founder_muts
        # method to avoid dataframes.
        # Make a dataframe of possible mutations to the parent, add a column giving the
        # seqeunce context of each mutation, and initialize columns to record counts
        # and branch lengths.
        counts_df, parent_seq = self.poss_muts.possible_muts_from_founder_muts(p_muts)
        parent_motif = lambda s: parent_seq[s - 2 : s + 1]
        counts_df["parent_motif"] = counts_df.site.apply(parent_motif)
        counts_df["actual_count"] = 0
        counts_df["branch_length"] = 0
        pcps = (parent.id, parent_seq, [])

        # Iterate over the children of the parent node, and count mutations to each
        n_filtered = 0
        for node in parent.children:
            
            # Get the mutations on the branch and only analyze them if they
            # pass the filter
            nt_muts = node.mutations
            codon_muts = self.translations[node.id]
            if self.fails_filters(nt_muts, codon_muts):
                continue
            n_filtered += 1
            pcps[2].append((node.id, nt_muts))
            counts_df["branch_length"] += len(nt_muts)
            to_increment = counts_df.query("nt_mut.isin(@node.mutations)").index
            counts_df.loc[to_increment, "actual_count"] += 1

        # If at least some branches pass the filter, then record data
        # TODO: ask Chris about logic
        if n_filtered != 0:
            records = counts_df.to_records(index=False).tolist()
            counts_dict = Counter({x[:-2]: x[-2] for x in records})
            branches_dict = Counter({x[:-2]: x[-1] for x in records})
        else:
            counts_dict, branches_dict = Counter(), Counter()

        return n_filtered, counts_dict, branches_dict, pcps

    def mut_counters_to_df(self, actual_counts, branch_lengths):
        """
        Return a DataFrame with data from two counters like those returned by
        self.count_mutations_from_parent.
        """
        zipped = zip(actual_counts.items(), branch_lengths.values())
        data = [(*k, v, w) for (k, v), w in zipped]
        columns = [*self.key_names, "actual_count", "branch_length"]
        the_df = pd.DataFrame(data, columns=columns)
        return the_df

    def site_counter_to_df(self, counter, n_filtered):
        """
        Return a dataframe with site conservation information from a Counter.
        """
        the_df = pd.DataFrame.from_dict(counter, orient="index")
        the_df.reset_index(inplace=True)
        the_df.rename(columns={"index": "site", 0: "mut_count"}, inplace=True)
        the_df["frac_mut"] = the_df.mut_count / n_filtered
        the_df["frac_cons"] = 1 - the_df.frac_mut
        the_df.sort_values("mut_count", ascending=False, inplace=True)
        return the_df


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

        # Iterate over all nodes
        for i, node in enumerate(self.int_nodes, 1):
            if i % 1000 == 0:
                print(f"processing node {i}; {t0+time():0.1f} seconds")

            # If the node is a parent node, count mutations going to each of its children
            # TODO: check with Chris about the line I added and the Counter logic
            # TODO: update n_filtered variable name
            if node.children:
                count, actual_counts, branch_lens, pcps = self.count_mutations_from_parent(node)
                n_filtered += count
                all_actual_counts.update(actual_counts)
                all_branch_lengths.update(branch_lens)
                if count != 0:
                    all_pcps.append(pcps)

        # Aggregate the counts and branch lengths across all nodes
        all_counts_df = self.mut_counters_to_df(all_actual_counts, all_branch_lengths)
        all_pcps_df = self.pcp_list_to_df(all_pcps)

        # Add metadata to the counts dataframe
        all_counts_df['mut_type']  = all_counts_df['wt_nt'] + all_counts_df['mut_nt']

        def get_mut_class(wt_aa, mut_aa):
            if wt_aa == mut_aa:
                return 'synonymous'
            elif mut_aa == '*':
                return 'nonsense'
            else:
                return 'nonsynonymous'
        all_counts_df['mut_class'] = all_counts_df.apply(lambda row: get_mut_class(row['wt_aa'], row['mut_aa']), axis=1)

        # Compress the counts dataframe to have one row per site
        all_counts_df['parent_motif'] = all_counts_df['parent_motif'].fillna('N.A.')
        groupby_cols = ['site', 'nt_mut', 'wt_nt', 'mut_nt', 'parent_motif', 'actual_count', 'branch_length']
        all_counts_df = all_counts_df.groupby(groupby_cols).agg({
            'gene': lambda x: ';'.join(str(val) for val in x),
            'codon_position': lambda x: ';'.join(str(val) for val in x),
            'codon_site': lambda x: ';'.join(str(val) for val in x),
            'wt_codon': lambda x: ';'.join(str(val) for val in x),
            'mut_codon': lambda x: ';'.join(str(val) for val in x),
            'wt_aa': lambda x: ';'.join(str(val) for val in x),
            'mut_aa': lambda x: ';'.join(str(val) for val in x),
            'aa_mut': lambda x: ';'.join(str(val) for val in x),
            'mut_class': lambda x: ';'.join(str(val) for val in x),
        }).reset_index()

        return n_filtered, all_counts_df, all_pcps_df

    # def parallel_count_mutations_on_tree(self, processes=8, batch_size=1000):
    #     """
    #     ...would be nice, but not there yet...
    #     """
    #     raise NotImplementedError(
    #         "Can't multiprocess with bte tree, need to do extra book keeping."
    #     )
    #     n_filtered, all_actual_counts, all_branch_lengths = 0, Counter(), Counter()
    #     batch_count = ceil(self.n_int_nodes / batch_size)
    #     t0 = -time()
    #     with Pool(processes=processes) as pool:
    #         for b in range(batch_count):
    #             if b == batch_count - 1:
    #                 nodes = self.int_nodes[b * batch_size : (b + 1) * batch_size]
    #             else:
    #                 nodes = self.int_nodes[b * batch_size :]
    #             batch_result = pool.map(self.count_mutations_from_parent, nodes)
    #             for count, actual_counts, branch_lens in batch_result:
    #                 n_filtered += count
    #                 all_actual_counts.update(actual_counts)
    #                 all_branch_lengths.update(branch_lens)

    #         print(
    #             f"{b} batches of {batch_size} nodes processed in {t0+time():0.1f} seconds"
    #         )
    #     all_counts_df = self.mut_counters_to_df(all_actual_counts, all_branch_lengths)
    #     return n_filtered, all_counts_df

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
        the_df = pd.DataFrame(rows, columns=cols)
        return the_df

    def make_and_write_all_counts(
        self,
        all_counts_path,
        all_pcps_path,
    ):
        """Calculate various counts and write to file."""

        print(f"Number of nodes in tree: {self.n_nodes}")
        print(f"Number of internal nodes in tree: {self.n_int_nodes}")

        print(f"Counting mutations on tree...")
        n_filtered, all_counts_df, all_pcps_df = self.count_mutations_on_tree()
        print(f"Number of nodes passing filter: {n_filtered}")

        all_counts_df.to_csv(all_counts_path, index=False)
        all_pcps_df.to_csv(all_pcps_path, index=False)

        return None

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Count mutations along a phylogenetic tree")
    parser.add_argument('--tree_path', required=True, help='Path to the tree file')
    parser.add_argument('--coding_site_path', required=True, help='Path to coding sites CSV file')
    parser.add_argument('--fasta_path', required=True, help='Path to reference FASTA file')
    parser.add_argument('--gtf_path', required=True, help='Path to GTF annotation file')
    parser.add_argument('--all_counts_path', required=True, help='Output path for mutation counts CSV')
    parser.add_argument('--all_pcps_path', required=True, help='Output path for parent-child pairs CSV')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Initialize an instance of CountsHelper
    counts_helper = CountsHelper(
        args.tree_path,
        args.coding_site_path,
        args.fasta_path,
        args.gtf_path
    )
    
    # Compute counts and write them to files
    counts_helper.make_and_write_all_counts(
        args.all_counts_path,
        args.all_pcps_path
    )

if __name__ == "__main__":
    main()