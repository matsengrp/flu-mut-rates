"""
Tests for the origin-sequence guards in ExpectedCalc.

Run from the repository root with `pytest scripts/`.
"""

import doctest

import pytest

import ExpectedCalc
from ExpectedCalc import apply_muts, check_tree_origin, parse_haplotype_muts

# A short stand-in for a tree's origin sequence. Site n is ORIGIN[n - 1].
ORIGIN = "ATGCATGCAT"

# Maps a nucleotide to one it is not, for building deliberately disagreeing mutations.
WRONG_NT = {"A": "C", "C": "G", "G": "T", "T": "A"}


class Node:
    """
    Minimal stand-in for a bte.MATNode: mutations on the branch into this node, in
    reflocalt format and relative to its parent.
    """

    def __init__(self, id, mutations=(), children=()):
        self.id = id
        self.mutations = list(mutations)
        self.children = list(children)


# --- parse_haplotype_muts: checks one node's accumulated haplotype ---


def test_parses_agreeing_haplotype():
    assert parse_haplotype_muts(["A1G", "G3T", "T10C"], ORIGIN) == {
        1: "G",
        3: "T",
        10: "C",
    }


def test_empty_haplotype_is_no_mutations():
    assert parse_haplotype_muts([], ORIGIN) == {}


def test_accepts_a_set_as_bte_returns_one():
    assert parse_haplotype_muts({"C4A", "A5T"}, ORIGIN) == {4: "A", 5: "T"}


def test_parsed_muts_reconstruct_the_sequence():
    muts = parse_haplotype_muts(["A1G", "T10C"], ORIGIN)
    assert apply_muts(ORIGIN, muts) == "GTGCATGCAC"


def test_disagreeing_parent_base_raises():
    # Site 3 is 'G' in ORIGIN, so a mutation claiming 'A' there is a different origin.
    with pytest.raises(ValueError, match=r"site 3 is 'G' in the origin but 'A'"):
        parse_haplotype_muts(["A1G", "A3T"], ORIGIN)


def test_error_reports_how_many_mutations_disagree():
    with pytest.raises(ValueError, match=r"2 of 3 haplotype mutations disagree"):
        parse_haplotype_muts(["A1G", "A3T", "A4G"], ORIGIN)


def test_error_truncates_long_mismatch_lists():
    haplotype = [f"{WRONG_NT[ORIGIN[i]]}{i + 1}G" for i in range(8)]
    with pytest.raises(ValueError, match=r"8 of 8 .*and 3 more"):
        parse_haplotype_muts(haplotype, ORIGIN)


def test_site_past_end_of_origin_raises():
    with pytest.raises(ValueError, match=r"outside the origin sequence of length 10"):
        parse_haplotype_muts(["A11G"], ORIGIN)


def test_site_zero_raises():
    # Sites are 1-indexed, so site 0 would silently read the last base of the origin.
    with pytest.raises(ValueError, match=r"outside the origin sequence"):
        parse_haplotype_muts(["A0G"], ORIGIN)


def test_wrong_origin_at_one_site_is_caught():
    """
    The failure mode this guard exists for: an origin fasta that differs from the tree's
    real origin at a single site. Every haplotype covering that site must be rejected.
    """
    wrong_origin = "C" + ORIGIN[1:]
    assert sum(a != b for a, b in zip(ORIGIN, wrong_origin)) == 1
    haplotype = ["A1G", "T2G"]
    parse_haplotype_muts(haplotype, ORIGIN)  # fine against the real origin
    with pytest.raises(ValueError, match=r"site 1 is 'C' in the origin but 'A'"):
        parse_haplotype_muts(haplotype, wrong_origin)


# --- check_tree_origin: checks every branch's mutations in one walk ---


def consistent_tree():
    """
    A tree whose mutations agree with ORIGIN ("ATGCATGCAT"). Site 1 mutates twice down
    the left lineage (A->G->C), so the second mutation's parent base is the intermediate
    'G' rather than the origin's 'A'.
    """
    return Node(
        "root",
        [],
        [
            Node("a", ["A1G", "T2C"], [Node("a1", ["G1C"]), Node("a2", ["G3A"])]),
            Node("b", ["C4T"], [Node("b1", ["A5C"])]),
        ],
    )


def test_consistent_tree_passes_and_counts_mutations():
    assert check_tree_origin(consistent_tree(), ORIGIN) == 6


def test_empty_tree_passes():
    assert check_tree_origin(Node("root"), ORIGIN) == 0


def test_sibling_lineages_do_not_see_each_others_state():
    """
    The walk must undo a branch's mutations when it backtracks. Node 'b' is a sibling of
    'a', so site 1 is still the origin's 'A' on b's branch even though 'a' changed it.
    """
    tree = Node(
        "root",
        [],
        [Node("a", ["A1G"], [Node("a1", ["G1T"])]), Node("b", ["A1C"])],
    )
    assert check_tree_origin(tree, ORIGIN) == 3


def test_stale_sibling_state_is_caught():
    # 'b' claims a parent base of 'G' at site 1, which only holds inside a's lineage.
    tree = Node("root", [], [Node("a", ["A1G"]), Node("b", ["G1C"])])
    with pytest.raises(ValueError, match=r"1 of 2 mutations in the tree disagree"):
        check_tree_origin(tree, ORIGIN)


def test_wrong_origin_is_caught_and_sites_reported():
    # Swap sites 1 and 3 of the origin: the tree's mutations there no longer line up.
    wrong = "C" + ORIGIN[1] + "A" + ORIGIN[3:]
    with pytest.raises(ValueError, match=r"at 2 site\(s\): 1, 3"):
        check_tree_origin(consistent_tree(), wrong)


def test_error_truncates_long_site_lists():
    # One child per site 1-8, each claiming a parent base the origin does not have.
    tree = Node(
        "root",
        [],
        [Node(f"n{i}", [f"{WRONG_NT[ORIGIN[i]]}{i + 1}A"]) for i in range(8)],
    )
    with pytest.raises(ValueError, match=r"at 8 site\(s\): 1, 2, 3, 4, 5, and 3 more"):
        check_tree_origin(tree, ORIGIN, max_reported=5)


def test_mutation_past_end_of_origin_raises():
    tree = Node("root", [], [Node("a", ["A11G"])])
    with pytest.raises(ValueError, match=r"Node 'a' has mutation 'A11G' at site 11"):
        check_tree_origin(tree, ORIGIN)


def test_parse_haplotype_muts_docstring_examples():
    tests = [
        t
        for t in doctest.DocTestFinder().find(ExpectedCalc.parse_haplotype_muts)
        if t.examples
    ]
    assert tests, "the guard's docstring should carry runnable examples"
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    for test in tests:
        runner.run(test)
    assert runner.failures == 0
