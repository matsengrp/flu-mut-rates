"""
Tests for the origin-sequence guards in ExpectedCalc.

Run from the repository root with `pytest scripts/`.
"""

import doctest
import sys

import pytest

import ExpectedCalc
from ExpectedCalc import (
    apply_muts,
    check_tree_origin,
    parse_haplotype_muts,
    parse_reflocalt,
    summarize_disagreements,
)

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
    with pytest.raises(
        ValueError, match=r"at site 0, outside the origin sequence of length 10"
    ):
        parse_haplotype_muts(["A0G"], ORIGIN)


def test_duplicate_site_in_haplotype_raises():
    """
    A haplotype is stated relative to one origin, so a site cannot appear twice. The old
    dict-comprehension parse silently kept whichever entry came last.
    """
    with pytest.raises(ValueError, match=r"lists site 3 more than once"):
        parse_haplotype_muts(["G3T", "G3A"], ORIGIN)


def test_malformed_mutation_in_haplotype_raises_a_clear_error():
    with pytest.raises(ValueError, match=r"Malformed mutation string 'G3'"):
        parse_haplotype_muts(["G3"], ORIGIN)


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


def test_same_site_twice_on_one_branch_is_consistent():
    """
    Regression: one branch may carry two mutations at the same site (A->G->C on the way
    into a single node). The undo entries must then be unwound newest-first, or the
    intermediate 'G' leaks out of the branch and a later sibling that legitimately
    mutates site 1 from the origin's 'A' is wrongly rejected.
    """
    tree = Node(
        "root",
        [],
        [Node("a", ["A1G", "G1C"]), Node("b", ["A1T"])],
    )
    assert check_tree_origin(tree, ORIGIN) == 3


def test_same_site_twice_on_one_branch_still_catches_a_real_break():
    # The middle step is wrong: after A1G site 1 is 'G', so a mutation from 'T' is bogus.
    tree = Node("root", [], [Node("a", ["A1G", "T1C"])])
    with pytest.raises(ValueError, match=r"1 of 2 mutations in the tree disagree"):
        check_tree_origin(tree, ORIGIN)


def test_stale_sibling_state_is_caught():
    # 'b' claims a parent base of 'G' at site 1, which only holds inside a's lineage.
    tree = Node("root", [], [Node("a", ["A1G"]), Node("b", ["G1C"])])
    with pytest.raises(ValueError, match=r"1 of 2 mutations in the tree disagree"):
        check_tree_origin(tree, ORIGIN)


def test_children_are_visited_in_list_order():
    """
    Children are pushed onto the LIFO stack in reverse so they are visited in
    `node.children` order. The first-seen offending node is the one reported, so this
    pins the order down rather than leaving it to stack mechanics.
    """
    tree = Node("root", [], [Node("first", ["G1A"]), Node("second", ["G1A"])])
    with pytest.raises(ValueError, match=r"e\.g\. node 'first'"):
        check_tree_origin(tree, ORIGIN)


def test_wrong_origin_is_caught_and_sites_reported():
    # Change sites 1 and 3 of the origin: the tree's mutations there no longer line up.
    wrong = "C" + ORIGIN[1] + "A" + ORIGIN[3:]
    with pytest.raises(
        ValueError, match=r"at 2 site\(s\): site 1 \(e\.g\. node 'a'\), site 3"
    ):
        check_tree_origin(consistent_tree(), wrong)


def test_error_truncates_long_site_lists():
    # One child per site 1-8, each claiming a parent base the origin does not have.
    tree = Node(
        "root",
        [],
        [Node(f"n{i}", [f"{WRONG_NT[ORIGIN[i]]}{i + 1}A"]) for i in range(8)],
    )
    with pytest.raises(ValueError, match=r"at 8 site\(s\): site 1 .*, and 3 more"):
        check_tree_origin(tree, ORIGIN, max_reported=5)


def test_no_more_suffix_at_the_truncation_boundary():
    # Exactly max_reported disagreeing sites: every one is named, nothing is elided.
    tree = Node(
        "root",
        [],
        [Node(f"n{i}", [f"{WRONG_NT[ORIGIN[i]]}{i + 1}A"]) for i in range(5)],
    )
    with pytest.raises(ValueError) as exc:
        check_tree_origin(tree, ORIGIN, max_reported=5)
    assert "more" not in str(exc.value).split("The origin sequence")[0]


def test_mutation_past_end_of_origin_raises():
    tree = Node("root", [], [Node("a", ["A11G"])])
    with pytest.raises(ValueError, match=r"Node 'a' has mutation 'A11G' at site 11"):
        check_tree_origin(tree, ORIGIN)


def test_malformed_mutation_in_tree_raises_a_clear_error():
    tree = Node("root", [], [Node("a", ["A1"])])
    with pytest.raises(ValueError, match=r"Malformed mutation string 'A1'"):
        check_tree_origin(tree, ORIGIN)


def test_deep_linear_chain_does_not_hit_a_recursion_limit():
    """
    Real flu MATs are deep and ladder-shaped, which is why the walk is iterative. A chain
    far longer than Python's recursion limit must still be checked without blowing up.
    """
    depth = sys.getrecursionlimit() * 3
    node = Node("leaf", ["A5C"])
    for i in range(depth):
        node = Node(f"n{i}", ["T2C"] if i % 2 else ["C2T"], [node])
    assert check_tree_origin(Node("root", [], [node]), ORIGIN) == depth + 1


# --- shared parsing / message helpers ---


def test_parse_reflocalt_splits_a_well_formed_string():
    assert parse_reflocalt("A123G") == (123, "A", "G")


@pytest.mark.parametrize("bad", ["", "A", "AG", "A1GG", "AA1G", "A12", "123", "A-1G"])
def test_parse_reflocalt_rejects_malformed_strings(bad):
    # 'A12' is the dangerous one: sliced blindly it yields a mutant "base" of '2'.
    with pytest.raises(ValueError, match=r"Malformed mutation string"):
        parse_reflocalt(bad)


def test_summarize_disagreements_elides_and_counts():
    assert summarize_disagreements([1, 2, 3, 4], 2) == "1, 2, and 2 more"
    assert summarize_disagreements([1, 2], 2) == "1, 2"
    assert summarize_disagreements([1], 5) == "1"
    assert summarize_disagreements(["a", "b", "c"], 2, sep="; ") == "a; b; and 1 more"


@pytest.mark.parametrize(
    "func",
    [
        ExpectedCalc.parse_haplotype_muts,
        ExpectedCalc.parse_reflocalt,
        ExpectedCalc.summarize_disagreements,
    ],
    ids=lambda f: f.__name__,
)
def test_docstring_examples(func):
    """
    `DocTestRunner.run` records failures rather than raising, so assert on the count.
    """
    tests = [t for t in doctest.DocTestFinder().find(func) if t.examples]
    assert tests, f"{func.__name__} should carry runnable docstring examples"
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    for test in tests:
        runner.run(test, out=lambda s: None)
    assert runner.failures == 0, f"{runner.failures} doctest failure(s)"
