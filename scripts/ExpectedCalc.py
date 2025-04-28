"""
A class and helpers for computing expected counts across founder subtrees.
"""

import pandas as pd
import numpy as np
import re
from collections import defaultdict


# This version will work with a table of rates.


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
    """
    Provided a path to a fasta file containing a single sequence, load that sequence
    and return it as a string.
    """
    with open(reference_file, "r") as fh:
        next(fh)
        reference_seq = "".join(map(str.strip, fh))
    return reference_seq


def all_poss_muts(site_map):
    """
    Given a site map dataframe with rows 'site', and 'wt_nt', return an exploded
    dataframe with a row for all possible mutations to the input nts.
    """
    ret = site_map.copy()
    alt_nuc = [
        {"A": "C", "T": "C", "C": "T", "G": "C"},
        {"A": "T", "T": "A", "C": "A", "G": "T"},
        {"A": "G", "T": "G", "C": "G", "G": "A"},
    ]
    re_pattern = re.compile("A|T|C|G")
    for j in range(3):
        fn = lambda x: alt_nuc[j][x.group(0)]
        ret[j + 1] = ret.wt_nt.str.replace(re_pattern, fn, regex=True)

    ret_tall = ret.melt(id_vars=["site", "wt_nt"], value_name="mut_nt").drop(
        "variable", axis=1
    )
    return ret_tall


def apply_muts(sequence, muts):
    """
    Return the sequence after applying mutations.

    Args:
        sequence (str): The sequence.
        muts (dict): A dictionary of mutations. The keys are integers for the site of a
            mutation (site indices start at 1); the values are strings for the resulting
            nucleotide at the site. For example, muts[25] = 'G' indicates a mutation
            resulting in 'G' at site 25 (which is at index 24 of sequence).
    """
    if len(muts) == 0:
        return sequence
    sites = sorted(muts.keys())
    substrings = [sequence[: sites[0] - 1]]
    for start, stop in zip(sites[:-1], sites[1:]):
        tobase = muts[start]
        seq = sequence[start : stop - 1]
        substrings.append(tobase)
        substrings.append(seq)
    substrings.append(muts[sites[-1]])
    substrings.append(sequence[sites[-1] :])

    return "".join(substrings)


def wide_codon_sites_to_tall(coding_sites):
    """
    Given a dataframe like PossibleMutations.coding_sites, where each row lists all
    genes for a site, return a new dataframe like PossibleMutations.tall_coding_sites,
    where sites for overlapping genes are are listed as multiple rows.
    """
    # Maybe not the most efficient way, but this function is only called once (at
    # instantiation of PossibleMutations).
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
    return pd.DataFrame(coding_sites_tall)


class PossibleMutations:
    """
    Given a reference sequence (fasta path) and coding sites table (csv path, see
    "ref_coding_sites.py"), this first computes all possible nucleotide mutations
    (triple the reference genome size), then uses the annotation data to translate
    codons such that for each of the possible nucleotide mutations, you know the
    resulting wildtype translation, and mutation translation.

    Attributes:
        bad_site_list (list): A list of low quality nucleotide sites that are omitted in
            all calculations.
        coding_sites (pd.DataFrame): A dataframe with one row per nucleotide site that
            lists: the codon position (1, 2, or 3); codon site (position in the
            gene); and gene for the nucleotide site. The values for sites in overlapping
            genes are listed in a single row as semi-colon separated values.
        gene_fix (bool): When True, remove duplicate entries in self.coding_sites. See
           issue 46 (https://github.com/matsengrp/s2trajectory/issues/46) for details.
        ref_seq (str): The unmutated reference genome.
        reference_mutations_df (pd.DataFrame): A dataframe describing all possible
            nucleotide site mutations (of the reference sequence) in terms of the codon
            information in self.tall_coding_sites. Bad sites are omitted...
        reference_site_map (pd.DataFrame): The reference sequence as a dataframe listing
            site and nucleotide. Bad sites are omitted...
        tall_coding_sites (pd.DataFrame): A dataframe with the same information as
            self.coding_sites. The values for sites in overlapping genes are listed in
            separate rows.


    Example
    -------
    Intialized with the SCV2 Wuhan reference strain and coding sites,
    we will have access to the ``reference_mutations_df`` attribute:

    >>> poss = PossibleMutations(fasta_path, coding_sites, secondary_struct, True)
    >>> poss.reference_mutations_df.query("site.isin([13467, 13468, 13469])")
            site wt_nt mut_nt   gene  codon_position codon_sites wt_codon mut_codon wt_aa mut_aa  is_synonymous      nuc      aa            ID
    40398  13467     A      C  ORF1a               2        4401      AAC       ACC     N      T          False  A13467C  N4401T  ORF1a-N4401T
    40399  13467     A      T  ORF1a               2        4401      AAC       ATC     N      I          False  A13467T  N4401I  ORF1a-N4401I
    40400  13467     A      G  ORF1a               2        4401      AAC       AGC     N      S          False  A13467G  N4401S  ORF1a-N4401S
    40401  13468     C      T  ORF1a               3        4401      AAC       AAT     N      N           True  C13468T  N4401N  ORF1a-N4401N
    40402  13468     C      T  ORF1b               1        4402      CGG       TGG     R      W          False  C13468T  R4402W  ORF1b-R4402W
    40403  13468     C      A  ORF1a               3        4401      AAC       AAA     N      K          False  C13468A  N4401K  ORF1a-N4401K
    40404  13468     C      A  ORF1b               1        4402      CGG       AGG     R      R           True  C13468A  R4402R  ORF1b-R4402R
    40405  13468     C      G  ORF1a               3        4401      AAC       AAG     N      K          False  C13468G  N4401K  ORF1a-N4401K
    40406  13468     C      G  ORF1b               1        4402      CGG       GGG     R      G          False  C13468G  R4402G  ORF1b-R4402G
    40407  13469     G      C  ORF1a               1        4402      GGG       CGG     G      R          False  G13469C  G4402R  ORF1a-G4402R
    40408  13469     G      C  ORF1b               2        4402      CGG       CCG     R      P          False  G13469C  R4402P  ORF1b-R4402P
    40409  13469     G      T  ORF1a               1        4402      GGG       TGG     G      W          False  G13469T  G4402W  ORF1a-G4402W
    40410  13469     G      T  ORF1b               2        4402      CGG       CTG     R      L          False  G13469T  R4402L  ORF1b-R4402L
    40411  13469     G      A  ORF1a               1        4402      GGG       AGG     G      R          False  G13469A  G4402R  ORF1a-G4402R
    40412  13469     G      A  ORF1b               2        4402      CGG       CAG     R      Q          False  G13469A  R4402Q  ORF1b-R4402Q

    Above, you can see the entries for just three reference sites (13467, 13468, 13469)
    in the reference sequence. Sites 13468 and 13469 are in the overlap of the ORF1a and
    ORF1b genes, and so we make an entry for each possible codon mutation in a gene.
    For example, the nucleotide mutation C->T at site 13468 in ORF1a is the synonymous
    amino acid mutation N->N, while in ORF1b it is the non-synonymous R->W.

    Additionally, the class provides the ability to obtain a modified version of the
    above table for new subtree founder variants. For example, we have the possible
    mutations after the nucleotide at site 13468 has mutated to an A.

    >>> founder_muts_df = poss.possible_muts_from_founder_muts({13468: "A"})
    >>> founder_muts_df.query("site.isin([13467, 13468, 13469])")
            site wt_nt mut_nt   gene  codon_position codon_sites wt_codon mut_codon wt_aa mut_aa  is_synonymous      nuc      aa            ID
    88584  13467     A      C  ORF1a               2        4401      AAA       ACA     K      T          False  A13467C  K4401T  ORF1a-K4401T
    88585  13467     A      T  ORF1a               2        4401      AAA       ATA     K      I          False  A13467T  K4401I  ORF1a-K4401I
    88586  13467     A      G  ORF1a               2        4401      AAA       AGA     K      R          False  A13467G  K4401R  ORF1a-K4401R
    88587  13468     A      C  ORF1a               3        4401      AAA       AAC     K      N          False  A13468C  K4401N  ORF1a-K4401N
    88588  13468     A      C  ORF1b               1        4402      AGG       CGG     R      R           True  A13468C  R4402R  ORF1b-R4402R
    88589  13468     A      T  ORF1a               3        4401      AAA       AAT     K      N          False  A13468T  K4401N  ORF1a-K4401N
    88590  13468     A      T  ORF1b               1        4402      AGG       TGG     R      W          False  A13468T  R4402W  ORF1b-R4402W
    88591  13468     A      G  ORF1a               3        4401      AAA       AAG     K      K           True  A13468G  K4401K  ORF1a-K4401K
    88592  13468     A      G  ORF1b               1        4402      AGG       GGG     R      G          False  A13468G  R4402G  ORF1b-R4402G
    88593  13469     G      C  ORF1a               1        4402      GGG       CGG     G      R          False  G13469C  G4402R  ORF1a-G4402R
    88594  13469     G      C  ORF1b               2        4402      AGG       ACG     R      T          False  G13469C  R4402T  ORF1b-R4402T
    88595  13469     G      T  ORF1a               1        4402      GGG       TGG     G      W          False  G13469T  G4402W  ORF1a-G4402W
    88596  13469     G      T  ORF1b               2        4402      AGG       ATG     R      M          False  G13469T  R4402M  ORF1b-R4402M
    88597  13469     G      A  ORF1a               1        4402      GGG       AGG     G      R          False  G13469A  G4402R  ORF1a-G4402R
    88598  13469     G      A  ORF1b               2        4402      AGG       AAG     R      K          False  G13469A  R4402K  ORF1b-R4402K
    """  # noqa: E501

    def __init__(
        self,
        ref_seq: str,
        coding_sites: str,
        bad_sites_path: str,
        secondary_struct: str,
        rates_table: str,
        gene_fix=True,
    ):
        self.gene_fix = gene_fix
        self.ref_seq = load_reference_from_fasta(ref_seq)
        self.basepair_df = pd.read_csv(secondary_struct)
        self.rates_df = pd.read_csv(rates_table)
        self.boundary_df = self.rates_df.groupby("mut_type").nt_site_boundary.first()
        # Now that we've recorded the site boundaries, we can drop them from rates_df.
        self.rates_df.drop(columns=["nt_site_boundary"], inplace=True)

        self.set_coding_sites_df(coding_sites)
        if len(self.coding_sites) != len(self.ref_seq):
            raise ValueError(
                f"reference sequence and codon sites must be the same length"
                f", got {len(self.ref_seq)}, and {len(self.coding_sites)}"
            )
        if not np.all(self.coding_sites.site == np.arange(1, len(self.ref_seq) + 1)):
            raise ValueError("codon sites must be in order of reference sequence")

        self.tall_coding_sites = wide_codon_sites_to_tall(self.coding_sites)
        self.load_bad_sites_from_file(bad_sites_path)
        self.reference_site_map = pd.DataFrame(
            {"site": np.arange(1, len(self.ref_seq) + 1), "wt_nt": list(self.ref_seq)}
        )
        self.reference_site_map.query("~site.isin(@self.bad_site_list)", inplace=True)
        self.reference_mutations_df = self.possible_mutations_df()

    def load_bad_sites_from_file(self, bad_sites_path):
        """Load a list of bad sites from file."""
        with open(bad_sites_path) as the_file:
            self.bad_site_list = sorted(map(int, map(str.strip, the_file)))
        return None

    def add_extra_featues(self, the_df, ref_seq):
        """
        Given a dataframe with columns "site", "wt_nt", and "mut_nt" and a reference
        sequence, return a new dataframe with additional columns  "mut_type", "motif"
        (3-mer local sequence context), "unpaired" (secondary rna structure),
        "nt_site_boundary" (the site boundary for global context), and
        "nt_site_before_boundary" (whether or not the site is before the boundary).
        """
        motif_at_site = lambda s: ref_seq[s - 2 : s + 1]
        the_df = the_df.merge(self.basepair_df, on="site")
        the_df["motif"] = the_df.site.apply(motif_at_site)
        the_df["mut_type"] = the_df.wt_nt + the_df.mut_nt
        the_df = the_df.merge(self.boundary_df, on="mut_type")
        the_df["nt_site_before_boundary"] = the_df.site < the_df.nt_site_boundary

        return the_df

    def set_coding_sites_df(self, file_path):
        """
        Current fix for issue 46 (https://github.com/matsengrp/s2trajectory/issues/46).
        Once that is fixed upstream, need only read in csv.
        """
        coding_sites = pd.read_csv(file_path, keep_default_na=False)
        if self.gene_fix:
            to_fix = coding_sites.query("266 <= site <=13467")
            codon_position = to_fix.codon_position.str.split(";").str.get(0)
            codon_site = to_fix.codon_site.str.split(";").str.get(0)
            coding_sites.loc[to_fix.index, "codon_position"] = codon_position
            coding_sites.loc[to_fix.index, "codon_site"] = codon_site
            coding_sites.loc[to_fix.index, "gene"] = "ORF1a"

            to_fix = coding_sites.query("site==13468")
            coding_sites.loc[to_fix.index, "codon_position"] = "3;1"
            coding_sites.loc[to_fix.index, "codon_site"] = "4401;4402"
            coding_sites.loc[to_fix.index, "gene"] = "ORF1a;ORF1b"

            to_fix = coding_sites.query("13469 <= site <= 13480")
            coding_sites.loc[to_fix.index, "gene"] = "ORF1a;ORF1b"

            to_fix = coding_sites.query("13481 <= site <= 21552")
            coding_sites.loc[to_fix.index, "gene"] = "ORF1b"

        self.coding_sites = coding_sites
        return None

    def possible_mutations_df(self, site_map=None, ref_seq=None):
        """
        Return a new dataframe fitting the description of self.reference_mutations_df,
        based on the given reference sequence and restricted to entries of site_map.

        Parameters:
            site_map (pd.DataFrame): A dataframe with columns "site" and "wt_nt". When
                None, default to self.reference_site_map.
            ref_seq (str): The reference sequence. When None, default to self.ref_seq.

        See attributes documentation for details.
        """
        if site_map is None or len(site_map) == 0:
            site_map = self.reference_site_map
        ref_seq = self.ref_seq if ref_seq is None else ref_seq

        query = f"site.isin({site_map.site.to_list()})"
        total_annotated_sites = self.tall_coding_sites.query(query)
        all_muts = all_poss_muts(site_map)
        all_muts = all_muts.merge(total_annotated_sites, on="site", how="outer")

        # Drop rows with noncoding codon positions now, as these rows give nans.
        query = "codon_position == 'noncoding'"
        all_muts.drop(index=all_muts.query(query).index, inplace=True)
        all_muts.codon_position = all_muts.codon_position.astype(int)

        # Is it possible to vectorize the wt_codon calculation?
        all_muts["wt_codon"] = (all_muts.site - all_muts.codon_position).apply(
            lambda x: ref_seq[x : x + 3]
        )
        mut_codons = []
        for j in range(1, 4):
            rows = all_muts.query("codon_position==@j")
            slice = rows.wt_codon.str.slice
            mut_codons.append(slice(stop=j - 1) + rows.mut_nt + slice(start=j))
        all_muts["mut_codon"] = pd.concat(mut_codons)
        all_muts["wt_aa"] = all_muts.wt_codon.map(CODON_TABLE)
        all_muts["mut_aa"] = all_muts.mut_codon.map(CODON_TABLE)
        all_muts["is_synonymous"] = all_muts.wt_aa == all_muts.mut_aa
        all_muts["nuc"] = all_muts.wt_nt + all_muts.site.astype(str) + all_muts.mut_nt
        all_muts["aa"] = (
            all_muts.wt_aa + all_muts.codon_sites.astype(str) + all_muts.mut_aa
        )
        all_muts["ID"] = all_muts.gene + "-" + all_muts.aa

        return all_muts

    def possible_muts_from_founder_muts(self, founder_muts):
        """
        Return a new dataframe constructed by applying mutations to
        self.reference_mutations_df.

        Parameters:
            founder_muts (dict): A dictionary of mutations. The keys are integers for
                the site of a mutation (site indices start at 1); the values are strings
                for the resulting nucleotide at the site. For example, muts[25] = 'G'
                is a mutation resulting in 'G' at site 25 (which is at index 24 of
                sequence).
        """
        # cols_to_keep = ["site", "wt_nt", "mut_nt", "is_synonymous", "ID"]
        # ret = self.reference_mutations_df[cols_to_keep].copy()
        ret = self.reference_mutations_df.copy()
        if len(founder_muts) == 0:
            return ret, self.ref_seq

        # Take all sites with nucleotide mutations, expand to all sites in codons with
        # nucleotide mutations, and determine nucleotides (after mutation) at the sites.
        mut_sites = list(founder_muts.keys())
        query = "site.isin(@mut_sites) & codon_position != 'noncoding'"
        codon_data = self.tall_coding_sites.query(query)
        start_positions = codon_data.site - codon_data.codon_position.astype(int) + 1
        sites = sorted({start + j for start in start_positions for j in range(3)})
        wt_nts = [
            founder_muts[site] if site in founder_muts else self.ref_seq[site - 1]
            for site in sites
        ]
        new_site_map = {"site": sites, "wt_nt": wt_nts}

        # Drop all rows with sites that are being redefined by the site map.
        to_drop = ret.query(f"site.isin({new_site_map['site']})").index
        ret.drop(to_drop, inplace=True)

        # Apply muts to get new reference_sequence.
        new_ref_seq = apply_muts(self.ref_seq, founder_muts)

        # get the new possible mutations rows for the new site map
        new_poss_muts = self.possible_mutations_df(
            pd.DataFrame(new_site_map), new_ref_seq
        )

        # Append those to the ret, and return.
        ret = pd.concat([ret, new_poss_muts], ignore_index=True)

        return ret, new_ref_seq

    def get_actual_expected_for_subtree(
        self, counts, report_nucleotides=False, context_aware=True
    ):
        """
        Given a dataframe of counts for a single subtree, compute the actual and
        expected counts for every possible mutation to the founder sequence. Return a
        dataframe with columns "actual_count" and "expected_count", indexed by the
        codon-mutation ID (<gene>-<wtaa><codonsite><mutaa>). Thread safe.

        Parameters:
            count (pd.DataFrame): The dataframe of actual mutation counts.
            context_aware (bool): When True, the expected counts are based on rates
                of synonmous mutations grouped by mutation type, motif, secondary
                rna structure, and global sequence position. These rates are computed
                in another project. When False, the expected counts are instead the
                average number of synonmyous mutations grouped by mutation type in the
                given subtree.
        """
        founder_mutations = counts.mutations.iat[0]
        # Get all possible mutations from the founder.
        all_possible_founder_muts, ref_seq = self.possible_muts_from_founder_muts(
            founder_mutations
        )
        all_possible_founder_muts["actual_count"] = 0
        # Stack observed counts with all possible
        all_possible_counts = pd.concat([all_possible_founder_muts, counts])

        cols_to_keep = {"actual_count"}
        if context_aware:
            # Add sequence context feature columns.
            all_possible_counts = self.add_extra_featues(all_possible_counts, ref_seq)
            group_keys = [
                "ID",
                "site",
                "mut_type",
                "wt_nt",
                "mut_nt",
                "motif",
                "unpaired",
                "nt_site_before_boundary",
                "is_synonymous",
            ]
        else:
            group_keys = ["ID", "site", "wt_nt", "mut_nt", "is_synonymous"]
        cols_to_keep.update(group_keys)
        cols_to_drop = set(all_possible_counts.columns) - cols_to_keep
        all_possible_counts.drop(columns=cols_to_drop, inplace=True)

        # Combine the actual observed counts with the 0 counts from all possible founder
        # mutations.
        all_possible_counts = all_possible_counts.groupby(
            group_keys, as_index=False
        ).sum()

        # Compute expected counts.
        if context_aware:
            subtree_mut_count = counts.actual_count.sum()
            merge_keys = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
            all_possible_counts = all_possible_counts.merge(
                self.rates_df, on=merge_keys
            )
            all_possible_counts.rename(columns={"rate": "expected_count"}, inplace=True)
            all_possible_counts.expected_count *= subtree_mut_count
        else:
            expected_counts = (
                all_possible_counts.query("is_synonymous")
                .groupby(by=["wt_nt", "mut_nt"])
                .agg(expected_count=("actual_count", "mean"))
                .reset_index()
            )
            all_possible_counts = all_possible_counts.merge(
                expected_counts, on=["wt_nt", "mut_nt"]
            )

        if report_nucleotides:
            nuc_counts = all_possible_counts.copy()
            nuc_counts["syn"] = "syn-"
            non_syn_entries = nuc_counts.query("~is_synonymous").index
            nuc_counts.loc[non_syn_entries, "syn"] = "nonsyn-"
            nuc_counts["ID"] = (
                "nuc-"
                + nuc_counts.syn
                + nuc_counts.wt_nt
                + nuc_counts.site.astype("str")
                + nuc_counts.mut_nt
            )
            all_possible_counts = pd.concat((all_possible_counts, nuc_counts))

        cols_to_keep = {"ID", "actual_count", "expected_count"}
        cols_to_drop = set(all_possible_counts.columns) - cols_to_keep
        all_possible_counts.drop(columns=cols_to_drop, inplace=True)

        all_possible_counts = all_possible_counts.groupby("ID").sum()

        return all_possible_counts
