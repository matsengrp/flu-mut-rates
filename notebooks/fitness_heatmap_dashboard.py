import marimo

__generated_with = "0.20.4"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import altair as alt
    import gffutils
    from Bio import SeqIO
    from Bio.Seq import Seq
    alt.data_transformers.disable_max_rows()
    return Seq, SeqIO, alt, gffutils, mo, np, pd


@app.cell
def _(pd):
    # Load and preprocess fitness effects data (computed once at startup)
    _raw = pd.read_csv("../results/aa_fitness_effects.csv", keep_default_na=False)

    # Keep only host=="all" mutations (all classes); drop PB1-F2 fusions
    _df = _raw[
        (_raw["host"] == "all") &
        (_raw["mut_class"].isin(["nonsynonymous", "synonymous", "nonsense"])) &
        (~_raw["gene"].str.contains("PB1-F2"))
    ].copy()

    # Normalize gene to uppercase
    _df["gene"] = _df["gene"].str.upper()

    # Explode overlapping-ORF rows (e.g. "NEP;NS1") on all four joined columns
    for _col in ["gene", "codon_site", "wt_aa", "mut_aa"]:
        _df[_col] = _df[_col].str.split(";")
    _df = _df.explode(["gene", "codon_site", "wt_aa", "mut_aa"])

    # After explosion, drop rows that are spuriously synonymous in this gene
    # (e.g. a nonsynonymous/nonsense mutation that happens to be synonymous in
    # the other overlapping ORF).  True synonymous mutations (mut_class ==
    # "synonymous") always have wt_aa == mut_aa, so we keep those.
    _df = _df[
        (_df["mut_class"] == "synonymous") | (_df["wt_aa"] != _df["mut_aa"])
    ].copy()
    _df["codon_site"] = _df["codon_site"].astype(int)

    preprocessed_df = _df.reset_index(drop=True)
    return (preprocessed_df,)


@app.cell
def _():
    import yaml
    with open("../config.yaml") as _f:
        _config = yaml.safe_load(_f)
    # TODO: Switch to results/{segment}/{subtype} once flu-usher pipeline is updated
    FLU_USHER_RESULTS = _config["data_dir"]

    GENE_TO_SEGMENT = {
        "HA": "HA", "NA": "NA", "NP": "NP",
        "PA": "PA", "PB1": "PB1", "PB2": "PB2",
        "M1": "MP", "M2": "MP",
        "NS1": "NS", "NEP": "NS",
    }
    return FLU_USHER_RESULTS, GENE_TO_SEGMENT


@app.cell
def _(FLU_USHER_RESULTS, Seq, SeqIO, gffutils):
    def load_reference_aa(segment, subtype, gene):
        """Return dict {codon_site (int): amino_acid (str)} for the given gene.

        Uses the same CDS concatenation logic as make_coding_sites.py:
        CDS features are sorted by start position, concatenated in order, then
        translated. codon_site is 1-indexed into the resulting protein.
        """
        base = f"{FLU_USHER_RESULTS}/{segment}/{subtype}"
        db = gffutils.create_db(
            f"{base}/curated_reference.gff",
            ":memory:",
            force=True,
            merge_strategy="create_unique",
        )
        ref_seq = str(next(SeqIO.parse(f"{base}/curated_reference.fasta", "fasta")).seq)

        # Collect CDS exons for the requested gene
        cdss = []
        for cds in db.features_of_type("CDS"):
            if cds.attributes.get("gene", [""])[0].upper() == gene.upper():
                cdss.append((cds.start, cds.end))
        cdss.sort()

        if not cdss:
            return {}

        # Concatenate exon sequences and translate (trim to multiple of 3 if needed)
        cds_nt = "".join(ref_seq[s - 1 : e] for s, e in cdss)
        cds_nt = cds_nt[: len(cds_nt) - len(cds_nt) % 3]
        protein = str(Seq(cds_nt).translate(to_stop=True))

        return {i + 1: aa for i, aa in enumerate(protein)}

    return (load_reference_aa,)


@app.cell
def _(mo):
    HA_SUBTYPES = ["H1", "H3", "H5", "H7", "H9"]
    NA_SUBTYPES = ["N1", "N2", "N6", "N8", "N9"]
    PROTEINS = ["HA", "NA", "NP", "PA", "PB1", "PB2", "M1", "M2", "NS1", "NEP"]

    protein = mo.ui.dropdown(options=PROTEINS, value="HA", label="Protein")
    return HA_SUBTYPES, NA_SUBTYPES, protein


@app.cell
def _(HA_SUBTYPES, NA_SUBTYPES, mo, protein):
    if protein.value == "HA":
        _opts, _default = HA_SUBTYPES, "H1"
    elif protein.value == "NA":
        _opts, _default = NA_SUBTYPES, "N1"
    else:
        _opts, _default = ["all"], "all"

    subtype = mo.ui.dropdown(options=_opts, value=_default, label="Subtype")
    return (subtype,)


@app.cell
def _(mo):
    min_expected_count = mo.ui.slider(
        start=0, stop=200, step=5, value=10,
        label="Min expected count", show_value=True,
    )
    min_actual_count = mo.ui.slider(
        start=0, stop=200, step=5, value=0,
        label="Min actual count", show_value=True,
    )
    count_filter_mode = mo.ui.radio(
        options=["either (expected OR actual)", "both (expected AND actual)"],
        value="either (expected OR actual)",
        label="Count filter mode",
    )
    return count_filter_mode, min_actual_count, min_expected_count


@app.cell
def _(GENE_TO_SEGMENT, load_reference_aa, protein, subtype):
    _segment = GENE_TO_SEGMENT.get(protein.value, protein.value)
    _ref_subtype = subtype.value if protein.value in ("HA", "NA") else "all"
    reference_aa = load_reference_aa(_segment, _ref_subtype, protein.value)
    return (reference_aa,)


@app.cell
def _(
    count_filter_mode,
    min_actual_count,
    min_expected_count,
    preprocessed_df,
    protein,
    reference_aa,
    subtype,
):
    _df = preprocessed_df[preprocessed_df["gene"] == protein.value].copy()

    # Subtype filter
    if protein.value in ("HA", "NA"):
        _df = _df[_df["subtype"] == subtype.value]
    else:
        _df = _df[_df["subtype"] == "all"]

    # Count filter
    _min_exp = min_expected_count.value
    _min_act = min_actual_count.value
    if "OR" in count_filter_mode.value:
        _count_mask = (_df["expected_count"] >= _min_exp) | (_df["actual_count"] >= _min_act)
    else:
        _count_mask = (_df["expected_count"] >= _min_exp) & (_df["actual_count"] >= _min_act)
    _df = _df[_count_mask]

    # Reference filter: only show mutations where wt_aa matches the reference
    # sequence amino acid (applied to all mutation classes)
    if reference_aa:
        _df = _df[_df["codon_site"].map(reference_aa) == _df["wt_aa"]]

    plot_data = _df.copy()
    return (plot_data,)


@app.cell
def _(alt, np, plot_data):
    AA_ORDER = [
        "K", "R", "H",             # positively charged
        "D", "E",                  # negatively charged
        "N", "Q", "S", "T",        # polar uncharged
        "P", "G", "A", "V", "L", "I", "M",  # nonpolar aliphatic
        "F", "W", "Y",             # aromatic
        "C",                       # cysteine
        "*",                       # stop
    ]

    color_scale = alt.Scale(
        scheme="redblue",
        domainMid=0,
        domainMin=-6,
        domainMax=2,
        clamp=True,
    )

    # Brush on quantitative x so numeric interval selection works correctly
    site_brush = alt.selection_interval(
        encodings=["x"],
        mark=alt.BrushConfig(stroke="black", strokeWidth=2),
    )

    # --- Panel 1: site zoom bar ---
    _sites_df = plot_data[["codon_site"]].drop_duplicates().sort_values("codon_site")
    site_zoom_bar = (
        alt.Chart(_sites_df)
        .mark_rect(color="steelblue", opacity=0.6)
        .encode(x=alt.X("codon_site:Q", title=None, axis=None))
        .add_params(site_brush)
        .properties(width=900, height=20, title="Drag to zoom")
    )

    # --- Panel 2: site-level line plot (mean nonsynonymous fitness per site) ---
    _site_agg = (
        plot_data[plot_data["mut_class"] == "nonsynonymous"]
        .groupby("codon_site", as_index=False)["delta_fitness"]
        .mean()
        .rename(columns={"delta_fitness": "mean_fitness"})
    )
    lineplot = (
        alt.Chart(_site_agg)
        .mark_line(point=True, opacity=0.8)
        .transform_filter(site_brush)
        .encode(
            x=alt.X("codon_site:Q", title="Site"),
            y=alt.Y("mean_fitness:Q", title="Mean nonsynonymous fitness effect",
                    scale=alt.Scale(zero=False)),
            tooltip=[
                alt.Tooltip("codon_site:Q", title="Site"),
                alt.Tooltip("mean_fitness:Q", title="Mean nonsynonymous fitness effect", format=".3f"),
            ],
        )
        .properties(width=900, height=130)
    )

    # --- Panel 3: per-mutation heatmap ---
    # Background: grey bar per site to indicate sites that have data
    _sites_only = plot_data[["codon_site"]].drop_duplicates()
    background = (
        alt.Chart(_sites_only)
        .mark_bar(color="lightgray", opacity=0.2, width={"band": 1.0})
        .transform_filter(site_brush)
        .encode(x=alt.X("codon_site:O", title="Site",
                         axis=alt.Axis(labelOverlap=True, labelAngle=0)))
    )

    # Colored rectangles — clip extreme fitness values for display
    _heatmap_data = plot_data.copy()
    _heatmap_data["delta_fitness"] = np.clip(_heatmap_data["delta_fitness"], -6, 2)
    heatmap_rects = (
        alt.Chart(_heatmap_data)
        .mark_rect()
        .transform_filter(site_brush)
        .encode(
            x=alt.X("codon_site:O", title="Site",
                    axis=alt.Axis(labelOverlap=True, labelAngle=0)),
            y=alt.Y("mut_aa:N", title="Mutant AA",
                    scale=alt.Scale(domain=AA_ORDER),
                    axis=alt.Axis(grid=True, tickBand="extent",
                                  gridColor="#333333", gridWidth=2)),
            color=alt.Color("delta_fitness:Q", scale=color_scale, title="Fitness effect"),
            tooltip=[
                alt.Tooltip("codon_site:O", title="Site"),
                alt.Tooltip("wt_aa:N", title="Ref AA"),
                alt.Tooltip("mut_aa:N", title="Mut AA"),
                alt.Tooltip("delta_fitness:Q", title="Fitness effect", format=".3f"),
                alt.Tooltip("actual_count:Q", title="Actual count", format=".1f"),
                alt.Tooltip("expected_count:Q", title="Expected count", format=".1f"),
            ],
        )
    )

    # Wildtype (reference) markers: "x" at each (site, ref_aa) position
    _wt_df = plot_data[["codon_site", "wt_aa"]].drop_duplicates()
    wt_marks = (
        alt.Chart(_wt_df)
        .mark_text(text="x", color="black", fontSize=7)
        .transform_filter(site_brush)
        .encode(
            x=alt.X("codon_site:O", axis=alt.Axis(labelOverlap=True, labelAngle=0)),
            y=alt.Y("wt_aa:N", scale=alt.Scale(domain=AA_ORDER)),
        )
    )

    heatmap = (background + heatmap_rects + wt_marks).properties(width=alt.Step(16), height=320)

    # --- Final composition ---
    chart = (
        alt.vconcat(site_zoom_bar, lineplot, heatmap, spacing=4)
        .resolve_scale(color="independent")
        .configure_axis(grid=False, labelOverlap="parity")
        .configure_view(stroke=None)
        .configure(padding=10)
    )
    return (chart,)


@app.cell
def _(
    chart,
    count_filter_mode,
    min_actual_count,
    min_expected_count,
    mo,
    plot_data,
    protein,
    subtype,
):
    _n_muts = len(plot_data)
    _n_sites = plot_data["codon_site"].nunique()
    summary = mo.md(f"**{_n_muts} mutations shown across {_n_sites} sites.**")

    if protein.value in ("HA", "NA"):
        controls = mo.hstack(
            [protein, subtype, min_expected_count, min_actual_count, count_filter_mode],
            gap=2,
        )
    else:
        controls = mo.hstack(
            [protein, min_expected_count, min_actual_count, count_filter_mode],
            gap=2,
        )

    mo.vstack([controls, summary, mo.ui.altair_chart(chart)])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
