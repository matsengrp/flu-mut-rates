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
        gtf_path,
        host_tsv_path=None,
        host_groups=None,
    ):
        """
        Parameters:
            tree_path (str): The tree.
            coding_site_path (str): The coding sites dataframe (nucleotide site to codon
                and gene).
            fasta_path (str): The fasta for the reference sequence.
            gtf_path (str):
            host_tsv_path (str): Optional path to a TSV with columns 'node' and
                'host_group' giving the host classification of each node (internal and
                leaf). When provided, only branches whose parent and child have matching
                unambiguous host classifications are counted, and the output gains a
                'host' column.
            host_groups (iterable[str] or None): Optional set of host labels to keep.
                Branches whose parent's host is not in this set are silently skipped, so
                the output contains rows only for the listed hosts.
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

        # Optional host annotation: drop ambiguous nodes (those appearing >1x) from
        # the lookup, but remember them so we can give a precise error if any
        # internal node turns out to lack an unambiguous host.
        if host_tsv_path:
            self.node_host, self.ambiguous_nodes = self._load_host_tsv(host_tsv_path)
        else:
            self.node_host = None
            self.ambiguous_nodes = set()

        self.host_groups = set(host_groups) if host_groups else None

    @staticmethod
    def _load_host_tsv(path):
        df = pd.read_csv(path, sep="\t")
        counts = df["node"].value_counts()
        ambiguous = set(counts[counts > 1].index)
        unambiguous = df[~df["node"].isin(ambiguous)]
        print(
            f"Loaded host TSV {path}: {len(unambiguous)} unambiguous nodes, "
            f"{len(ambiguous)} ambiguous nodes dropped."
        )
        return dict(zip(unambiguous["node"], unambiguous["host_group"])), ambiguous

    def ref_motif(self, s):
        """Return the 3-mer motif centered at the site. Site indices begin at one."""
        m = self.ref_seq[s - 2 : s + 1]
        return m if len(m) == 3 else pd.NA

    def fails_filters(self, nt_muts, codon_muts, max_mutations=4):
        """
        Returns a tuple (fails, reason) indicating whether mutations at a node fail the filter.
        The filter requires that there are between one and max_mutations mutations and no two
        mutations target the same codon.

        Parameters:
            nt_mut (list): List of nucleotide mutations, string entries like "A1234G".
            codon_muts (list): List of codon mutations, entries are bte.AAChange objects.
            max_mutations (int): Maximum number of mutations allowed (default: 4).

        Returns:
            tuple: (bool, str | None) - (True, reason) if fails, (False, None) if passes.
                   Reasons: "too_many_mutations", "zero_mutations", "duplicate_codons"
        """
        num_muts = len(nt_muts)

        if num_muts > max_mutations:
            return True, "too_many_mutations"
        if num_muts == 0:
            return True, "zero_mutations"

        gene_codon_pairs = {(mut.gene, mut.aa_index) for mut in codon_muts}
        gene_codon_nuc_triples = {(mut.gene, mut.aa_index, mut.nuc) for mut in codon_muts}
        if len(gene_codon_pairs) < len(gene_codon_nuc_triples):
            return True, "duplicate_codons"

        return False, None

    def count_mutations_from_parent(self, parent, max_mutations=4):
        """
        Count mutations from the parent node to all of its child nodes. Return four
        values, which are the number of child nodes that pass the filter;
        a Counter mapping a tuple of values for self.key_names to the number of times
        the associating mutation type occured; a Counter mapping a tuple of values for
        self.key_names to the sum of branch lengths where the mutation type is possible;
        a list of parent child pairs, where each entry is a triple with of the
        parent node id, the parent node sequence, and a list of pairs of child node id
        and child node sequence; and a dictionary of filter statistics.

        Parameters:
            parent: Parent node to process.
            max_mutations (int): Maximum number of mutations allowed (default: 4).
        """

        # If host stratification is active, every parent must have an unambiguous
        # host classification in the TSV. Parents whose host is unambiguous but not
        # in the configured host_groups are silently skipped.
        parent_host = None
        if self.node_host is not None:
            parent_host = self.node_host.get(parent.id)
            if parent_host is None:
                if parent.id in self.ambiguous_nodes:
                    raise ValueError(
                        f"Parent node {parent.id!r} has ambiguous host annotation "
                        f"(appears multiple times in the host TSV)."
                    )
                raise ValueError(
                    f"Parent node {parent.id!r} is missing from the host TSV."
                )
            if self.host_groups is not None and parent_host not in self.host_groups:
                return 0, pd.DataFrame(), (parent.id, "", []), {
                    'total_branches': 0,
                    'passing_branches': 0,
                    'total_mutations': 0,
                    'passing_mutations': 0,
                    'filtered_by_too_many': 0,
                    'filtered_by_zero': 0,
                    'filtered_by_duplicates': 0,
                }

        # Get the set of mutations in the partent node relative to the reference
        parent_node_haplotype = self.tree.get_haplotype(parent.id)
        p_muts = {int(mut[1:-1]): mut[-1] for mut in parent_node_haplotype}

        # Make a dataframe of possible mutations to the parent, add a column giving the
        # seqeunce context of each mutation, and initialize columns to record counts
        # and branch lengths.
        counts_df, parent_seq = self.poss_muts.possible_muts_from_founder_muts(p_muts)

        # Compress the counts dataframe to have one row per site, where values for sites that
        # had more than one row are concatenated with a semicolon.
        groupby_cols = [
            'site', 'nt_mut', 'wt_nt', 'mut_nt'
        ]
        counts_df = counts_df.groupby(groupby_cols).agg({
            'gene': lambda x: ';'.join(str(val) for val in x),
            'codon_position': lambda x: ';'.join(str(val) for val in x),
            'codon_site': lambda x: ';'.join(str(val) for val in x),
            'wt_codon': lambda x: ';'.join(str(val) for val in x),
            'mut_codon': lambda x: ';'.join(str(val) for val in x),
            'wt_aa': lambda x: ';'.join(str(val) for val in x),
            'mut_aa': lambda x: ';'.join(str(val) for val in x),
            'aa_mut': lambda x: ';'.join(str(val) for val in x),
        }).reset_index()

        # Add the parent motif and initialize the actual count and branch length columns
        def parent_motif(s):
            m = parent_seq[s - 2 : s + 1]
            return m if len(m) == 3 else pd.NA
        counts_df["parent_motif"] = counts_df.site.apply(parent_motif)
        counts_df["actual_count"] = 0
        pcps = (parent.id, parent_seq, [])
        n_passing_muts = 0

        # Initialize filter statistics
        filter_stats = {
            'total_branches': 0,
            'passing_branches': 0,
            'total_mutations': 0,
            'passing_mutations': 0,
            'filtered_by_too_many': 0,
            'filtered_by_zero': 0,
            'filtered_by_duplicates': 0
        }

        # Iterate over the children of the parent node, and count mutations on
        # the branch going to each child, if the branch passes the filters
        n_passing_filters = 0
        for node in parent.children:
            # When host stratification is active, only count branches whose child
            # has the same unambiguous host as the parent.
            if self.node_host is not None:
                child_host = self.node_host.get(node.id)
                if child_host != parent_host:
                    continue

            nt_muts = node.mutations
            codon_muts = self.translations[node.id]
            num_muts = len(nt_muts)

            # Update total statistics
            filter_stats['total_branches'] += 1
            filter_stats['total_mutations'] += num_muts

            # Check if branch passes filters
            fails, reason = self.fails_filters(nt_muts, codon_muts, max_mutations)
            if fails:
                # Update filtered statistics by reason
                if reason == "too_many_mutations":
                    filter_stats['filtered_by_too_many'] += num_muts
                elif reason == "zero_mutations":
                    filter_stats['filtered_by_zero'] += num_muts
                elif reason == "duplicate_codons":
                    filter_stats['filtered_by_duplicates'] += num_muts
                continue

            # Branch passes filters
            n_passing_filters += 1
            filter_stats['passing_branches'] += 1
            filter_stats['passing_mutations'] += num_muts

            pcps[2].append((node.id, nt_muts))
            n_passing_muts += len(nt_muts)
            to_increment = counts_df.query("nt_mut.isin(@node.mutations)").index
            counts_df.loc[to_increment, "actual_count"] += 1

        branch_length = counts_df['actual_count'].sum()
        assert branch_length == n_passing_muts, \
            f"branch_length ({branch_length}) != n_passing_muts ({n_passing_muts})"
        counts_df["branch_length"] = branch_length
        
        # Count synonymous mutations by checking wt_aa == mut_aa. This is valid because
        # we already filtered out branches where two mutations target the same codon, so
        # each nucleotide mutation affects a distinct codon and can be classified
        # independently as synonymous or non-synonymous.
        counts_df["syn_branch_length"] = counts_df[counts_df['wt_aa'] == counts_df['mut_aa']]['actual_count'].sum()

        if self.node_host is not None:
            counts_df["host"] = parent_host

        return n_passing_filters, counts_df, pcps, filter_stats

    def count_mutations_on_tree(self, max_mutations=4):
        """
        Count the number of mutations along the tree. Return the number of nodes passing
        the filter, a dataframe of per site nucleotide mutations, a dataframe of
        parent child pairs, and aggregate filter statistics.

        Parameters:
            max_mutations (int): Maximum number of mutations allowed (default: 4).
        """
        # Initialize variables
        all_counts_df = pd.DataFrame()
        n_passing_filters = 0
        all_pcps = []
        t0 = -time()

        # Initialize aggregate filter statistics
        aggregate_stats = {
            'total_branches': 0,
            'passing_branches': 0,
            'total_mutations': 0,
            'passing_mutations': 0,
            'filtered_by_too_many': 0,
            'filtered_by_zero': 0,
            'filtered_by_duplicates': 0
        }

        # Iterate over all internal nodes. For each, record counts along branches to children,
        # only considering branches that pass filters. Also record each PCP passing filters.
        for i, node in enumerate(self.int_nodes, 1):

            if i % 1000 == 0:
                print(f"processing internal node {i}; {t0+time():0.1f} seconds")
            n_passing_filters_i, counts_df, pcps, filter_stats = self.count_mutations_from_parent(node, max_mutations)
            n_passing_filters += n_passing_filters_i

            # Accumulate filter statistics
            for key in aggregate_stats:
                aggregate_stats[key] += filter_stats[key]

            # Skip the internal node if no branches to its children passed the filters
            # (or, when host_groups is set, if the parent's host was filtered out)
            if counts_df.empty or counts_df['branch_length'].sum() == 0:
                continue

            # Add PCPs to the big list for the whole tree
            all_pcps.append(pcps)

            # Add counts and branch lengths to the big dataframe for the whole tree
            if all_counts_df.empty:
                all_counts_df = counts_df
            else:
                if n_passing_filters == 0:
                    continue
                cols = [
                    'site', 'nt_mut', 'wt_nt', 'mut_nt',
                    'gene', 'codon_position', 'codon_site',
                    'wt_codon', 'mut_codon', 'wt_aa', 'mut_aa',
                    'aa_mut', 'parent_motif'
                ]
                if self.node_host is not None:
                    cols.append('host')
                all_counts_df = (
                    pd.concat([all_counts_df, counts_df])
                    .groupby(cols, as_index=False)
                    .agg({
                        'actual_count' : 'sum',
                        'branch_length' : 'sum',
                        'syn_branch_length' : 'sum'
                    })
                )

        # Add metadata to the counts dataframe
        def get_mut_class(wt_aa_list, mut_aa_list):
            
            # Split the amino acid lists and determine mutation class
            wt_aa_list = wt_aa_list.split(';')
            mut_aa_list = mut_aa_list.split(';')
            mut_class_list = []
            for wt_aa, mut_aa in zip(wt_aa_list, mut_aa_list):
                if wt_aa == mut_aa:
                    mut_class_list.append('synonymous')
                elif mut_aa == '*':
                    mut_class_list.append('nonsense')
                else:
                    mut_class_list.append('nonsynonymous')
            
            # Determine the overall mutation class
            if 'nonsense' in mut_class_list:
                return 'nonsense'
            elif 'nonsynonymous' in mut_class_list:
                return 'nonsynonymous'
            else:
                return 'synonymous'

        all_counts_df['mut_class'] = all_counts_df.apply(lambda row: get_mut_class(row['wt_aa'], row['mut_aa']), axis=1)
        all_counts_df['mut_type']  = all_counts_df['wt_nt'] + all_counts_df['mut_nt']

        # Aggregate the counts and branch lengths across all nodes
        all_pcps_df = self.pcp_list_to_df(all_pcps)

        return n_passing_filters, all_counts_df, all_pcps_df, aggregate_stats

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
        max_mutations=4,
    ):
        """
        Calculate various counts and write to file.

        Parameters:
            all_counts_path (str): Output path for mutation counts CSV.
            all_pcps_path (str): Output path for parent-child pairs CSV (optional).
            max_mutations (int): Maximum number of mutations allowed (default: 4).
        """

        print(f"Number of nodes in tree: {self.n_nodes}")
        print(f"Number of internal nodes in tree: {self.n_int_nodes}")

        print(f"Counting mutations on tree...")
        n_passing_filters, all_counts_df, all_pcps_df, filter_stats = self.count_mutations_on_tree(max_mutations)
        print(f"Number of nodes passing filter: {n_passing_filters}")

        # Print detailed filtering statistics
        print("\n" + "=" * 60)
        print("MUTATION FILTERING STATISTICS")
        print("=" * 60)
        print(f"Filter threshold: max_mutations = {max_mutations}")
        print()

        total_branches = filter_stats['total_branches']
        passing_branches = filter_stats['passing_branches']
        filtered_branches = total_branches - passing_branches

        total_mutations = filter_stats['total_mutations']
        passing_mutations = filter_stats['passing_mutations']
        filtered_mutations = total_mutations - passing_mutations

        print(f"Total branches examined:              {total_branches:,}")
        print(f"Branches passing filters:             {passing_branches:,} ({100*passing_branches/total_branches:.1f}%)" if total_branches > 0 else "Branches passing filters:             0 (0.0%)")
        print(f"Branches filtered:                    {filtered_branches:,} ({100*filtered_branches/total_branches:.1f}%)" if total_branches > 0 else "Branches filtered:                    0 (0.0%)")
        print()

        print(f"Total mutations (all branches):       {total_mutations:,}")
        print(f"Mutations on passing branches:        {passing_mutations:,} ({100*passing_mutations/total_mutations:.1f}%)" if total_mutations > 0 else "Mutations on passing branches:        0 (0.0%)")
        print(f"Mutations on filtered branches:       {filtered_mutations:,} ({100*filtered_mutations/total_mutations:.1f}%)" if total_mutations > 0 else "Mutations on filtered branches:       0 (0.0%)")
        print()

        print("Mutations filtered by reason:")
        print(f"  - Too many mutations (>{max_mutations}):        {filter_stats['filtered_by_too_many']:,} ({100*filter_stats['filtered_by_too_many']/total_mutations:.1f}%)" if total_mutations > 0 else f"  - Too many mutations (>{max_mutations}):        0 (0.0%)")
        print(f"  - Zero mutations:                   {filter_stats['filtered_by_zero']:,} ({100*filter_stats['filtered_by_zero']/total_mutations:.1f}%)" if total_mutations > 0 else "  - Zero mutations:                   0 (0.0%)")
        print(f"  - Duplicate codon targets:          {filter_stats['filtered_by_duplicates']:,} ({100*filter_stats['filtered_by_duplicates']/total_mutations:.1f}%)" if total_mutations > 0 else "  - Duplicate codon targets:          0 (0.0%)")
        print()

        # Sanity check
        sum_filtered = (filter_stats['filtered_by_too_many'] +
                       filter_stats['filtered_by_zero'] +
                       filter_stats['filtered_by_duplicates'])
        if sum_filtered != filtered_mutations:
            print(f"WARNING: Filtered mutation counts don't match! Sum of categories: {sum_filtered:,}, Expected: {filtered_mutations:,}")
        else:
            print("Sanity check passed: Sum of filtered categories matches total filtered mutations.")

        print("=" * 60 + "\n")

        all_counts_df.to_csv(all_counts_path, index=False)
        if all_pcps_path is not None:
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
    parser.add_argument('--all_pcps_path', required=False, default=None, help='Output path for parent-child pairs CSV (optional)')
    parser.add_argument('--host_tsv', required=False, default=None, help='Optional TSV mapping node IDs to host_group; enables host-stratified counts')
    parser.add_argument('--host_groups', required=False, nargs='+', default=None, help='Optional set of host labels to keep; only branches whose parent has one of these hosts are recorded')

    # Parse arguments
    args = parser.parse_args()

    # Initialize an instance of CountsHelper
    counts_helper = CountsHelper(
        args.tree_path,
        args.coding_site_path,
        args.fasta_path,
        args.gtf_path,
        host_tsv_path=args.host_tsv,
        host_groups=args.host_groups,
    )
    
    # Compute counts and write them to files
    counts_helper.make_and_write_all_counts(
        args.all_counts_path,
        args.all_pcps_path
    )

if __name__ == "__main__":
    main()