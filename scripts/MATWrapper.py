"""
This script breaks up aggregation of subtree mutations into two steps.

* The first step is the initialization of a MATWrapper object, which loads the protobuf tree and time information,
    does translation of mutations on the tree and computes sequences at every node in the tree. All of this takes awhile,
    but if you want to experiment with different window and step sizes, it need all be done only once.
* The second step is in running the `MATWrapper.aggregate_data` method, which iterates through all subtrees in all windows,
    and aggregates mutation information for each.

See the ```if __name__ == "__main__"``` block at the bottom of this script for an example of how to use, or just change
values as needed there.
"""
import bte

from functools import cache
from warnings import warn
import csv
from datetime import datetime
from collections import Counter, defaultdict


def get_time_dict(chronumental_tsv_file): # , invalid_strains):
    """returns a dictionary mapping node ids to estimated times.
    to get a node's time, look up the node's id in this dictionary.

    times are parsed as datetime objects, which is convenient for
    comparison and describing time by interval offsets"""
    data = {}
    with open(chronumental_tsv_file) as fh:
        rd = csv.reader(fh, delimiter="\t", quotechar='"')
        # skip first line of csv
        next(rd)
        for name, time in rd:
            # if name not in invalid_strains:
            data[name] = datetime.fromisoformat(time)

    return data


def traverse_mat(root_node, leaf_condition=lambda n: n.is_leaf(), mask=lambda n: False):
    """Implements a preorder traversal, with a custom stopping condition.

    The function leaf_condition should take a node and return true iff the traversal should treat that node as a leaf.

    The function mask should take a node and return true whenever a node should not be visited. This is an alternative way
    to stop the traversal on a condition.

    The traversal will begin at the provided root_node, wherever
    that node may reside in the ambient tree.

    Since the usher tree is shallow, this shouldn't run into recursion depth limitations.
    """
    yield root_node
    if not leaf_condition(root_node):
        for cnode in root_node.children:
            if not mask(cnode):
                yield from traverse_mat(cnode, leaf_condition=leaf_condition, mask=mask)


def iter_leaves(root_node, leaf_condition=lambda n: n.is_leaf(), mask=lambda n: False):
    """Visits leaves as in the `traverse_mat` function, where a node can only be
    considered a leaf if `leaf_condition` is True, and no nodes will be visited
    for which `mask(node)` is True.

    Perhaps this function should be called `iter_filtered_leaves`, because even when a node
    has no children for which `mask` returns False, it still won't be yielded unless
    `leaf_condition` is True on that node. This differs from the behavior of `traverse_mat`,
    which will yield such a node, even if `leaf_condition` returns False on it."""
    if leaf_condition(root_node):
        yield root_node
    else:
        for cnode in root_node.children:
            if not mask(cnode):
                yield from iter_leaves(cnode, leaf_condition=leaf_condition, mask=mask)


def load_reference_from_fasta(reference_file):
    """Provided a path to a fasta file containing a single sequence, load
    that sequence and return it as a string"""
    with open(reference_file, "r") as fh:
        next(fh)
        reference_seq = ""
        for line in fh:
            reference_seq += line.strip()
    return reference_seq


def apply_muts(sequence, muts):
    """Apply the (one-based) mutations in muts (an iterable containing mutations
    in format like `A25G`) to the provided sequence, returning the result"""
    for mut in muts:
        frombase = mut[0]
        tobase = mut[-1]
        site = int(mut[1:-1]) - 1
        if sequence[site] != frombase:
            warn("parent base doesn't match existing base at the sequence")
        sequence = sequence[:site] + tobase + sequence[site + 1 :]
    return sequence


def get_codon(sequence, start_codon, end_codon, site):
    """For a 1-indexed site in sequence, and the start and end codons containing a mutation at that site,
    find the corresponding codon in the sequence"""
    offset = [i for i, (a, b) in enumerate(zip(start_codon, end_codon)) if a != b][0]
    return sequence[site - offset - 1 : site - offset + 2]


class MATWrapper:
    """This is a wrapper for a mutation annotated tree object plus timing information.
    On initialization, all data required to aggregate over time-window-defined subtrees
    is gathered (this can take some time), so that data can be aggregated for a given window and step size by
    calling the method `aggregate_data`."""

    def __init__(
        self,
        usher_tree_file,
        chronumental_dates_file,
        GTF_file,
        reference_fasta,
        # invalid_strains,
    ):
        self.tree = bte.MATree(usher_tree_file)
        print("loaded tree")
        self.nodes = self.tree.depth_first_expansion()
        #invalid_strains = []
        #with open(invalid_strains) as fh:
        #    for line in fh:
        #        invalid_strains.append(line.strip())
        #self.nodes = [
        #    node
        #    for node in self.tree.depth_first_expansion
        #    if node.id not in set(invalid_strains)
        #]
        self.ids_to_nodes = {node.id: node for node in self.nodes}
        self.times = get_time_dict(chronumental_dates_file)

        self.start_time = self.times[self.nodes[0].id]
        self.end_time = max(self.times.values())
        self.reference_seq = load_reference_from_fasta(reference_fasta)

        print("translating mutations")
        # nodes without mutations aren't represented in this dictionary, and even some
        # nodes with mutations don't appear to be represented either... not sure why
        # that is -- does translate not identify synonymous mutations? if that's the
        # case, how do we count those?
        #
        # this is a mapping from node id strings to lists of bte.AAChange records
        translations = self.tree.translate(GTF_file, reference_fasta)
        # assume nodes not in translations dictionary just don't have parent mutations:
        self.translations = defaultdict(lambda: [], translations)

    @cache
    def get_haplotype(self, node_id):
        return self.tree.get_haplotype(node_id)

    def export_branch_data(self, output_file):
        """Outputs a csv with columns `branch_index`, `start_time`,
        `end_time` and `mutation_list`. This includes data for all edges in the MAT,
        even those which are filtered by `aggregate_data`. Output is written to the filepath
        `output_file` in csv format."""

        rows = [("branch_index", "start_time", "end_time", "mutation_list")]
        for idx, node in enumerate(self.tree.depth_first_expansion()[1:]):
            rows.append(
                (
                    idx,
                    self.times[node.parent.id],
                    self.times[node.id],
                    " ".join(str(mut) for mut in node.mutations),
                )
            )

        with open(output_file, "w") as fh:
            writer = csv.writer(fh, delimiter=",")
            for row in rows:
                writer.writerow(row)

    def aggregate_data(self, window_size, step_size):
        """Given a window size and a step size, produce a list containing for each window start time
        a pair (start_time, subtrees_data), where subtrees_data is a list containing for each subtree
        in that window the MAT node ID and mutations in the founder sequence relative to the reference, and a list
        summarizing all mutations in that subtree (this list can be empty if the subtree contains
        no mutations).

        Each mutation is summarized as a tuple containing (in this order):
            0) gene name
            1) nucleotide index (1-based)
            2) amino acid index
            3) aa mutation (e.g. F274F)
            4) nucleotide mutation (e.g. C29095T)
            5) original codon (e.g. TCT)
            6) alternative codon (e.g. CCT)
            7) a boolean describing whether the mutation is synonymous
            8) a count of how many mutations in the subtree shared exactly the same data above ^^

        This completes Hugh's suggested implementation up to and including point 4.


        This method also sets self.subtree_records as a similarly-formatted list of lists of nodes belonging to
        each subtree in each window, with an ordering which matches the ordering of the data structure returned
        by this method. This makes it easy to examine individual subtrees, if desired.


        This data can also be output as two csv files using the `convert_subtree_data_to_dataframes` function or
        by providing the `-c` option to the cli. One csv will contain subtree founders and their mutations
        relative to the reference, and the other will contain a row for each mutation, containing the data listed above,
        in addition to window start time and founder node id.
        """

        def qc_filter(node, founder_seq, qc_counter=Counter()):
            """
            Return True if the branch above node should be considered, otherwise False.

            For use in eliminating branches that...
            1. Have more than four nucleotide mutations
            2. Contain more than one nucleotide mutation that is a reversion to the Wuhan-Hu-1
            reference sequence
            3. Contain more than one nucleotide mutation that was a reversion to the subtree
            founder sequence (is this necessary for us?)
            4. Contain more than one nucleotide mutation to the same codon
            """
            nuc_muts = node.mutations
            # 1:
            if len(nuc_muts) > 4:
                # raise ValueError(
                #    "OUI! node has more than four nucleotide mutations! This should have been filtered above"
                # )
                qc_counter.update({1: 1})
                return False
            # 2:
            if (
                sum(
                    1
                    for mut in nuc_muts
                    if mut[-1] == self.reference_seq[int(mut[1:-1]) - 1]
                )
                > 1
            ):
                qc_counter.update({2: 1})
                return False
            # 3:
            if (
                sum(1 for mut in nuc_muts if mut[-1] == founder_seq[int(mut[1:-1]) - 1])
                > 1
            ):
                qc_counter.update({3: 1})
                return False
            # 4:
            aa_muts = self.translations[node.id]
            if len(aa_muts) < 2:
                pass
            elif (
                max(Counter((aamut.gene, aamut.aa_index) for aamut in aa_muts).values())
                > 1
            ):
                qc_counter.update({4: 1})
                return False

            qc_counter.update({0: 1})
            return True

        # n_windows = (self.end_time - self.start_time - window_size) // step_size
        n_windows = (self.end_time - self.start_time) // step_size
        print("starting traversal")
        print(
            f"time range {self.start_time} to {self.end_time} with window size of {window_size} ({n_windows} steps)"
        )

        # identify founder nodes: for a time window t_i to t_i + window_size, a node is
        # a founder if it lies in the window, but its parent does not.
        # founders will be a list of pairs (start_time, founder_list), where
        # founder_list is the list of founder nodes in the window starting at
        # start_time.
        founders = []
        for time in [self.start_time + i * step_size for i in range(n_windows)]:

            def leaf_func(node):
                return self.times[node.id] >= time and (
                    node.parent is None or self.times[node.parent.id] < time
                )

            def mask_func(node):
                return self.times[node.id] > time + step_size

            founders.append(
                (time, list(iter_leaves(self.tree.root, leaf_func, mask_func)))
            )

        self.founders = founders
        # # TODO: building the list of founders as above is inefficient,
        # # because it requires beginning a tree traversal from the root for each
        # # window. If windows are non-overlapping we can do something like the
        # # following, where we complete only one full traversal of the tree:
        # for node_id, time in self.times.items():
        #     node = self.ids_to_nodes[node_id]
        #     window_idx = (time - self.start_time) // window_size
        #     window_start, founder_list = founders[window_idx]
        #     if node.parent is None or self.times[node.parent.id] < window_start:
        #             founder_list.append(node)

        # subtree_records will have (window start time, subtree_list) pairs, where the
        # subtree_list contains for each subtree in the window, a list of nodes in that
        # subtree (including the founder node as the first node in the list)
        print("populating subtrees")
        subtree_records = []
        for window_start, founder_list in founders:
            subtrees = []
            for founder_node in founder_list:
                nl = list(
                    traverse_mat(
                        founder_node,
                        mask=lambda n: self.times[n.id] > window_start + window_size,
                    )
                )
                # Check that subtree contains at least one node:
                if len(nl) > 1:
                    # Check that subtree contains at least one mutation, below the founder node:
                    # Founder node is first in nl node list
                    if sum(len(n.mutations) for n in nl[1:]) > 0:
                        subtrees.append(nl)
            subtree_records.append((window_start, subtrees))
        self.subtree_records = subtree_records

        qc_counter = Counter()

        def summarize_subtree(node_list):
            """Aggregate data about a single subtree whose nodes are contained in node_list, with
            the founder node the first in the list"""
            founder = node_list[0]
            internal_nodes = node_list[1:]
            subtree_size = len(node_list)
            branches_with_mutations = sum(
                1 for node in internal_nodes if len(node.mutations) > 0
            )
            founder_muts = self.get_haplotype(founder.id)
            founder_seq = apply_muts(self.reference_seq, founder_muts)
            records = []
            for node in filter(
                lambda n: qc_filter(n, founder_seq, qc_counter=qc_counter),
                internal_nodes,
            ):
                for aamut in self.translations[node.id]:
                    idx = aamut.nt_index
                    founder_codon = get_codon(
                        founder_seq,
                        aamut.original_codon,
                        aamut.alternative_codon,
                        aamut.nt_index,
                    )
                    if founder_codon == aamut.original_codon:
                        records.append(
                            (
                                aamut.gene,
                                aamut.nt_index,
                                aamut.aa_index,
                                aamut.aa,
                                aamut.nuc,
                                aamut.original_codon,
                                aamut.alternative_codon,
                                aamut.is_synonymous(),
                            )
                        )
            # This will deduplicate and add a 'count' column at the end
            records = [record + (count,) for record, count in Counter(records).items()]
            return (
                (founder.id, founder_muts, subtree_size, branches_with_mutations),
                records,
            )

        print("building subtree data")
        subtree_data = [
            (start_time, [summarize_subtree(node_list) for node_list in subtrees])
            for start_time, subtrees in subtree_records
        ]
        # To see how many branches were rejected for all the various reasons:
        # print(qc_counter)
        return subtree_data


def convert_subtree_data_to_dataframes(
    subtree_data, haplotypes_file, mutation_records_file
):
    """Write subtree mutation data to two csv files, named `haplotypes_file` and `mutation_records_file`.
    `haplotypes_file` will have five columns: `founder_id`, `window_start`, `subtree_size`, `branches_with_mutations`, and `mutations`,
    where `founder_id` is the MAT node ID of a subtree founder node, `subtree_size` is the number of nodes in the subtree,
    including the founder node, `branches_with_mutations` is the number of branches in the subtree with at least one mutation,
    and `mutations` is a string containing all mutations in that node's sequence, relative to the reference sequence
    (separated by spaces)."""
    founder_haplotypes = [
        (
            "founder_id",
            "window_start",
            "subtree_size",
            "branches_with_mutations",
            "mutations",
        )
    ]
    founder_haplotypes.extend(
        (
            founder_id,
            start_time,
            subtree_size,
            branches_with_mutations,
            " ".join(founder_muts),
        )
        for start_time, subtree_list in subtree_data
        for (
            founder_id,
            founder_muts,
            subtree_size,
            branches_with_mutations,
        ), records in subtree_list
    )
    mutation_records = [
        (
            "window_start",
            "founder_id",
            "gene",
            "nt_index",
            "aa_index",
            "aa",
            "nuc",
            "original_codon",
            "alternative_codon",
            "is_synonymous",
            "count",
        )
    ]
    mutation_records.extend(
        (str(start_time), founder_id) + record
        for start_time, subtree_list in subtree_data
        for (
            founder_id,
            founder_muts,
            subtree_size,
            branches_with_mutations,
        ), records in subtree_list
        for record in records
    )

    with open(haplotypes_file, "w") as fh:
        writer = csv.writer(fh, delimiter=",")
        for row in founder_haplotypes:
            writer.writerow(row)
    with open(mutation_records_file, "w") as fh:
        writer = csv.writer(fh, delimiter=",")
        for row in mutation_records:
            writer.writerow(row)
