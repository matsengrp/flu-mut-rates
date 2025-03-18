"""
A class and helpers for computing expected counts across founder subtrees.
"""

import pandas as pd
import numpy as np
from functools import lru_cache
from warnings import warn

from collections import defaultdict

# opt 1
CODON_TABLE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}


def load_reference_from_fasta(reference_file):
    """Provided a path to a fasta file containing a single sequence, load
    that sequence and return it as a string"""
    with open(reference_file, "r") as fh:
        next(fh)
        reference_seq = ""
        for line in fh:
            reference_seq += line.strip()
    return reference_seq


def all_poss_muts(site_map):  # , coding_sites_tall):
    """
    Given a site map dataframe with rows 'site', and 'wt_nt',
    return an exploded dataframe with a row
    for all possible mutations to the input nts.
    """
    # assert np.all(site_map.site == sorted(coding_sites_tall.site.unique()))

    nt_set = set(["A", "T", "C", "G"])
    ret = site_map.copy()
    # TODO padarallel this
    ret[[1, 2, 3]] = ret.apply(
        lambda x: tuple(nt_set.difference(set([x.wt_nt]))), axis=1, result_type="expand"
    )
    ret_tall = ret.melt(id_vars=["site", "wt_nt"], value_name="mut_nt").drop(
        "variable", axis=1
    )
    return ret_tall


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


def wide_codon_sites_to_tall(coding_sites):
    """
    convert a dataframe with columns that maps sites to each possible
    codon site, position, and gene, to a tall dataframe
    with columns 'site', 'codon_site', 'codon_pos', and 'gene'
    where sites are no longer unique.
    """
    coding_sites_tall = defaultdict(list)
    for row_idx, data in coding_sites.iterrows():
        genes = data.gene.split(";")
        codon_poss = data.codon_position.split(";")
        codon_sites = data.codon_site.split(";")
        for features in zip(genes, codon_poss, codon_sites):
            coding_sites_tall["site"].append(data.site)
            coding_sites_tall["gene"].append(features[0])
            coding_sites_tall["codon_position"].append(features[1])
            coding_sites_tall["codon_sites"].append(features[2])
    return coding_sites_tall


class PossibleMutations:
    """
    Given a reference sequence (fasta path) and coding sites table (csv path)
    (see ``ref_coding_sites.py``), this first computes all
    possible nucleotide mutations (3X reference genome size),
    then uses the annotation data to translate codons
    such that for each of the possible nt mutations,
    you know the resulting wildtype translation, and mutation translation.

    Example
    -------
    Intialized with the SCV2 Wuhan reference strain and coding sites,
    we will have access to the ``reference_mutations_df`` attribute:

    >>> poss = PossibleMutations(fasta_path, coding_sites)
    >>> poss.reference_mutations_df.query("site.isin([284, 285, 286])")
        site wt_nt mut_nt    gene codon_position codon_sites wt_codon mut_codon wt_aa mut_aa is_synonymous    nuc   aa          ID
    903   284     G      C   ORF1a              1           7      GGT       CGT     G      R         False  G284C  G7R   ORF1a-G7R
    904   284     G      C  ORF1ab              1           7      GGT       CGT     G      R         False  G284C  G7R  ORF1ab-G7R
    905   284     G      A   ORF1a              1           7      GGT       AGT     G      S         False  G284A  G7S   ORF1a-G7S
    906   284     G      A  ORF1ab              1           7      GGT       AGT     G      S         False  G284A  G7S  ORF1ab-G7S
    907   284     G      T   ORF1a              1           7      GGT       TGT     G      C         False  G284T  G7C   ORF1a-G7C
    908   284     G      T  ORF1ab              1           7      GGT       TGT     G      C         False  G284T  G7C  ORF1ab-G7C
    909   285     G      C   ORF1a              2           7      GGT       GCT     G      A         False  G285C  G7A   ORF1a-G7A
    910   285     G      C  ORF1ab              2           7      GGT       GCT     G      A         False  G285C  G7A  ORF1ab-G7A
    911   285     G      A   ORF1a              2           7      GGT       GAT     G      D         False  G285A  G7D   ORF1a-G7D
    912   285     G      A  ORF1ab              2           7      GGT       GAT     G      D         False  G285A  G7D  ORF1ab-G7D
    913   285     G      T   ORF1a              2           7      GGT       GTT     G      V         False  G285T  G7V   ORF1a-G7V
    914   285     G      T  ORF1ab              2           7      GGT       GTT     G      V         False  G285T  G7V  ORF1ab-G7V
    915   286     T      C   ORF1a              3           7      GGT       GGC     G      G          True  T286C  G7G   ORF1a-G7G
    916   286     T      C  ORF1ab              3           7      GGT       GGC     G      G          True  T286C  G7G  ORF1ab-G7G
    917   286     T      A   ORF1a              3           7      GGT       GGA     G      G          True  T286A  G7G   ORF1a-G7G
    918   286     T      A  ORF1ab              3           7      GGT       GGA     G      G          True  T286A  G7G  ORF1ab-G7G
    919   286     T      G   ORF1a              3           7      GGT       GGG     G      G          True  T286G  G7G   ORF1a-G7G
    920   286     T      G  ORF1ab              3           7      GGT       GGG     G      G          True  T286G  G7G  ORF1ab-G7G

    Above, you can see the entries for just three reference sites (284, 285, 286)
    in the reference sequence. These sites happen to exist within two overlapping
    codons for the ORF1a and ORF1ab genes, and so we make an entry for each possible
    _codon mutation_.

    Additionally, this class provides the ability to obtain a modified version
    of the reference mutations table above, for new subtree founder variants.

    >>> founder_muts_df = poss.possible_muts_from_founder_muts(tuple(["G285A"]))
    >>> founder_muts_df.query("site.isin([284, 285, 286])")
            site wt_nt mut_nt    gene codon_position codon_sites wt_codon mut_codon wt_aa mut_aa is_synonymous    nuc   aa          ID
    130215   284     G      C   ORF1a              1           7      GAT       CAT     D      H         False  G284C  D7H   ORF1a-D7H
    130216   284     G      C  ORF1ab              1           7      GAT       CAT     D      H         False  G284C  D7H  ORF1ab-D7H
    130217   284     G      A   ORF1a              1           7      GAT       AAT     D      N         False  G284A  D7N   ORF1a-D7N
    130218   284     G      A  ORF1ab              1           7      GAT       AAT     D      N         False  G284A  D7N  ORF1ab-D7N
    130219   284     G      T   ORF1a              1           7      GAT       TAT     D      Y         False  G284T  D7Y   ORF1a-D7Y
    130220   284     G      T  ORF1ab              1           7      GAT       TAT     D      Y         False  G284T  D7Y  ORF1ab-D7Y
    130221   285     A      C   ORF1a              2           7      GAT       GCT     D      A         False  A285C  D7A   ORF1a-D7A
    130222   285     A      C  ORF1ab              2           7      GAT       GCT     D      A         False  A285C  D7A  ORF1ab-D7A
    130223   285     A      T   ORF1a              2           7      GAT       GTT     D      V         False  A285T  D7V   ORF1a-D7V
    130224   285     A      T  ORF1ab              2           7      GAT       GTT     D      V         False  A285T  D7V  ORF1ab-D7V
    130225   285     A      G   ORF1a              2           7      GAT       GGT     D      G         False  A285G  D7G   ORF1a-D7G
    130226   285     A      G  ORF1ab              2           7      GAT       GGT     D      G         False  A285G  D7G  ORF1ab-D7G
    130227   286     T      C   ORF1a              3           7      GAT       GAC     D      D          True  T286C  D7D   ORF1a-D7D
    130228   286     T      C  ORF1ab              3           7      GAT       GAC     D      D          True  T286C  D7D  ORF1ab-D7D
    130229   286     T      A   ORF1a              3           7      GAT       GAA     D      E         False  T286A  D7E   ORF1a-D7E
    130230   286     T      A  ORF1ab              3           7      GAT       GAA     D      E         False  T286A  D7E  ORF1ab-D7E
    130231   286     T      G   ORF1a              3           7      GAT       GAG     D      E         False  T286G  D7E   ORF1a-D7E
    130232   286     T      G  ORF1ab              3           7      GAT       GAG     D      E         False  T286G  D7E  ORF1ab-D7E
    """  # noqa: E501

    def __init__(self, ref_seq: str, coding_sites: str):
        self.ref_seq = load_reference_from_fasta(ref_seq)
        self.coding_sites = pd.read_csv(coding_sites)

        if len(self.coding_sites) != len(self.ref_seq):
            raise ValueError(
                f"reference sequence and codon sites must be the same length"
                f", got {len(self.ref_seq)}, and {len(self.coding_sites)}"
            )

        if not np.all(self.coding_sites.site == np.arange(1, len(self.ref_seq) + 1)):
            raise ValueError("codon sites must be in order of reference sequence")
        self.tall_coding_sites = pd.DataFrame(
            wide_codon_sites_to_tall(self.coding_sites)
        )

        # get all possible mutations at each site
        self.reference_site_map = pd.DataFrame(
            {"site": np.arange(1, len(self.ref_seq) + 1), "wt_nt": list(self.ref_seq)}
        )

        self.reference_mutations_df = self.possible_mutations_df()

    def possible_mutations_df(self, site_map=None, ref_seq=None):
        """
        See class description.
        """

        site_map = self.reference_site_map if site_map is None else site_map
        ref_seq = self.ref_seq if ref_seq is None else ref_seq

        total_annotated_sites = self.tall_coding_sites.query(
            f"site.isin({site_map.site.to_list()})"
        )
        all_poss_muts_from_seq = all_poss_muts(site_map)

        ret_poss_muts = all_poss_muts_from_seq.merge(
            total_annotated_sites, on="site", how="outer"
        )

        def get_codon(site, codon_position, mut_nt, ref_seq):
            try:
                begin = int(site) - int(codon_position)
                wt_codon = ref_seq[begin : begin + 3]
                wt_codon_list = list(wt_codon)
                wt_codon_list[int(codon_position) - 1] = mut_nt
                mut_codon = "".join(wt_codon_list)
                # wt_aa, mut_aa = translate(wt_codon), translate(mut_codon)
                wt_aa, mut_aa = CODON_TABLE[wt_codon], CODON_TABLE[mut_codon]
                return wt_codon, mut_codon, wt_aa, mut_aa, wt_aa == mut_aa
            except ValueError:
                return np.nan, np.nan, np.nan, np.nan, np.nan

        # TODO pandarallel
        ret_poss_muts[
            ["wt_codon", "mut_codon", "wt_aa", "mut_aa", "is_synonymous"]
        ] = ret_poss_muts.apply(
            lambda x: get_codon(x.site, x.codon_position, x.mut_nt, ref_seq),
            axis=1,
            result_type="expand",
        )

        # TODO could just do this in the above apply, in parallel
        return (
            ret_poss_muts.assign(
                nuc=ret_poss_muts.wt_nt
                + ret_poss_muts.site.astype(str)
                + ret_poss_muts.mut_nt,
                aa=ret_poss_muts.wt_aa
                + ret_poss_muts.codon_sites.astype(str)
                + ret_poss_muts.mut_aa,
            )
            .assign(ID=lambda x: x.gene + "-" + x.aa)
            .dropna()
        )

    @lru_cache(maxsize=500)
    def possible_muts_from_founder_muts(self, founder_muts, slim=False):
        """
        Given a list of mutations in the format `A25G`, return a copy of the
        the reference_mutations_df dataframe with only the rows that correspond
        to those mutations.
        """
        # TODO dumb way to do this, but it works
        if founder_muts == tuple([]):
            ret = self.reference_mutations_df.copy()
            if slim:
                ret.drop(
                    [
                        "site",
                        # 'wt_nt',
                        # 'mut_nt',
                        "gene",
                        "codon_position",
                        "codon_sites",
                        "wt_codon",
                        "mut_codon",
                        "wt_aa",
                        "mut_aa",
                        # 'is_synonymous',
                        # 'nuc',
                        "aa",
                        # 'ID'
                    ],
                    inplace=True,
                    axis=1,
                )
                # ret.dropna(inplace=True)
            return ret

        new_site_map = defaultdict(list)
        new_site_map = {"site": [], "wt_nt": []}

        for mut in founder_muts:
            _wt, site, mutant = mut[0], int(mut[1:-1]), mut[-1]
            for idx, coding_sites in self.tall_coding_sites.query(
                f"site == {site}"
            ).iterrows():
                gene = coding_sites.gene  # noqa: F841
                codon_site = coding_sites.codon_sites  # noqa: F841
                codon_position = coding_sites.codon_position

                # TODO, could be faster by looking for ID's
                # instead of double conditional
                reference_muts_codon_df = self.reference_mutations_df.query(
                    "gene == @gene & codon_sites == @codon_site"
                )

                for codon_poss, poss_df in reference_muts_codon_df.groupby(
                    "codon_position"
                ):
                    nt_site = poss_df.site.values[0]
                    is_nt_mutation = codon_poss == codon_position
                    if nt_site in new_site_map["site"] and not is_nt_mutation:
                        # don't re-mutate a site that isn't our founder mutation
                        # if it's already been added, it should still be correct
                        continue
                    elif nt_site in new_site_map["site"] and is_nt_mutation:
                        # update sitemap
                        new_site_map["wt_nt"][
                            new_site_map["site"].index(nt_site)
                        ] = mutant
                    else:
                        new_site_map["site"].append(poss_df.site.values[0])
                        new_site_map["wt_nt"].append(
                            poss_df.wt_nt.values[0] if not is_nt_mutation else mutant
                        )
        print(new_site_map)

        ret = self.reference_mutations_df.copy()

        # drop all rows with sites that are being redined from the site map
        # reference_sites = ret.query(f"site.isin({sites})")
        # TODO do this inline
        to_drop = self.reference_mutations_df.query(
            f"site.isin({new_site_map['site']})"
        ).index

        ret.drop(to_drop, inplace=True)

        # apply muts to get new reference_sequence
        new_ref_seq = apply_muts(self.ref_seq, founder_muts)

        # get the new possible mutations rows for the new site map
        new_poss_muts = self.possible_mutations_df(
            pd.DataFrame(new_site_map), new_ref_seq
        )

        # append those to the ret, and return
        ret = pd.concat([ret, new_poss_muts], ignore_index=True)

        if slim:
            ret.drop(
                [
                    "site",
                    # 'wt_nt',
                    # 'mut_nt',
                    "gene",
                    "codon_position",
                    "codon_sites",
                    "wt_codon",
                    "mut_codon",
                    "wt_aa",
                    "mut_aa",
                    # 'is_synonymous',
                    # 'nuc',
                    "aa",
                    # 'ID'
                ],
                inplace=True,
                axis=1,
            )
            # ret.dropna(inplace=True)
        return ret


def get_actual_expected_from_subtree(
    counts, ecc, base_changes, window_founder_muts, founder_child_threshold
):
    """
    Given the a single group of subtree counts, and the respective
    founder mutations dataframe (must contain subtree founder entry),
    as well as an ExpectedCountsCalculator initialized from the reference
    sequence, compute the expected counts for every possible
    mutation to the founder sequence. Return a dataframe with columns
    Actual, and Expected, indexed by the codon-mutation ID
    (<gene>-<wtaa><codonsite><mutaa>).

    This function was designed to be parallelized across
    all trees in a given window.
    """

    # get unique founder
    uniq_founders = counts.founder_id.unique()
    assert len(uniq_founders) == 1
    founder = uniq_founders[0]
    subtree_founder_muts = window_founder_muts.query("founder_id == @founder")

    # threshold
    if subtree_founder_muts.subtree_size.values[0] <= founder_child_threshold:
        return pd.DataFrame(columns=["ID", "Actual", "Expected"]).astype(
            {"Actual": int, "Expected": float}
        )

    # get all possible mutations from a founder, assign count of zero
    all_possible_founder_muts = (
        ecc.possible_muts_from_founder_muts(
            tuple(subtree_founder_muts.mutations.values[0].split()), slim=True
        )
        .dropna()
        .assign(count=0, founder_id=founder)
    )
    all_possible_founder_muts["is_synonymous"] = all_possible_founder_muts[
        "is_synonymous"
    ].astype(bool)

    # merge the counts into all possible
    group_keys = ["ID", "nuc", "wt_nt", "mut_nt", "founder_id", "is_synonymous"]
    all_possible_counts = (
        pd.concat([all_possible_founder_muts, counts])
        .groupby(group_keys, as_index=False)
        .sum()
    )

    # opt 2
    # Make a dictionary of:
    #  key: (wt_nt, mut_nt)
    #  value:
    #    (
    #      total number of possible mutations with this base change,
    #      expected value for mut with this base change
    #    )
    E = {}
    for base_change in base_changes:
        # all possible mutations for this base change
        # Here we are doing an intermediary calculation to get the
        # total number of base changes for adding expected values to
        # the dataframe more quickly
        base_change_muts = all_possible_counts.query(
            f"wt_nt == '{base_change[0]}' & mut_nt == '{base_change[1]}'"
        )

        # now we can subset that to just synonymous mutations for
        # the expected calulation.
        synonymous_base_change_muts = base_change_muts.query("is_synonymous")

        E[base_change] = (
            base_change_muts["count"].shape[0],
            synonymous_base_change_muts["count"].mean(),
        )

    # now, since we know the total number of each mutation type, we can
    # assign values to the 'all possible' dataframe more quickly by sorting
    # rows in the df by the keys in the E dictionary, and then assigning the expected
    # values simply by repeating the expected value for each mutation type
    # into a flattened array.
    return (
        all_possible_counts.sort_values(["wt_nt", "mut_nt"])
        .assign(
            Expected=np.hstack([[E[key][1]] * E[key][0] for key in sorted(E.keys())])
        )
        .drop(
            [
                c
                for c in all_possible_counts.columns
                if c not in ["ID", "count", "Expected"]
            ],
            axis=1,
        )
        .rename({"count": "Actual"}, axis=1)
        .groupby("ID", as_index=False)
        .sum()
    )
