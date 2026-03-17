import marimo

__generated_with = "0.20.4"
app = marimo.App(width="full")


@app.cell
def _():
    import sys
    import warnings
    import marimo as mo
    import pandas as pd

    import altair as alt
    # Use default (inline JSON) transformer — avoids pyarrow dependency in WASM/Pyodide
    alt.data_transformers.enable("default")
    alt.data_transformers.disable_max_rows()
    # Suppress narwhals/altair compatibility warning (marimo wraps DataFrames in narwhals)
    warnings.filterwarnings("ignore", message="You passed a.*narwhals")
    # In WASM (Pyodide) fetch files over HTTP relative to the page URL;
    # locally they're one level up from the notebooks/ directory
    if "pyodide" in sys.modules:
        import js
        # Pyodide runs in a Web Worker; js.location.href is the worker script URL
        # (e.g. https://host/path/assets/worker.js). Strip from /assets/ onward
        # to recover the page base URL.
        _base = str(js.location.href).rsplit("/assets/", 1)[0]
        DATA_DIR = _base.rstrip("/") + "/results"
    else:
        DATA_DIR = "../results"
    return DATA_DIR, alt, mo, pd


@app.cell
def _(DATA_DIR, pd):
    import io, urllib.request as _urllib
    # Read via urllib + StringIO to prevent pandas double-decompressing gzip responses
    _url = f"{DATA_DIR}/nt_fitness_effects.csv"
    if _url.startswith("http"):
        with _urllib.urlopen(_url) as _r:
            df = pd.read_csv(io.StringIO(_r.read().decode("utf-8")), keep_default_na=False)
    else:
        df = pd.read_csv(_url, keep_default_na=False)
    return (df,)


@app.cell
def _(mo):
    HA_SUBTYPES = ["H1", "H3", "H5", "H7", "H9"]
    NA_SUBTYPES = ["N1", "N2", "N6", "N8", "N9"]
    SEGMENTS = ["HA", "NA", "NP", "PA", "PB1", "PB2", "MP", "NS"]
    segment = mo.ui.dropdown(options=SEGMENTS, value="HA", label="Segment")
    return HA_SUBTYPES, NA_SUBTYPES, segment


@app.cell
def _(HA_SUBTYPES, NA_SUBTYPES, mo, segment):
    if segment.value == "HA":
        _opts, _default = HA_SUBTYPES, "H1"
    elif segment.value == "NA":
        _opts, _default = NA_SUBTYPES, "N1"
    else:
        _opts, _default = ["all"], "all"
    subtype = mo.ui.dropdown(options=_opts, value=_default, label="Subtype")
    return (subtype,)


@app.cell
def _(mo):
    HOSTS = ["all", "human", "avian"]
    host = mo.ui.dropdown(options=HOSTS, value="all", label="Host")
    MUT_CLASSES = ["all", "nonsynonymous", "synonymous"]
    mut_class = mo.ui.dropdown(options=MUT_CLASSES, value="synonymous", label="Mutation class")
    min_count = mo.ui.slider(
        start=0, stop=200, step=5, value=10,
        label="Min count (actual or expected)", show_value=True,
    )
    return host, min_count, mut_class


@app.cell
def _(DATA_DIR):
    import json
    import urllib.request
    _url = f"{DATA_DIR}/reference_nt.json"
    if _url.startswith("http"):
        with urllib.request.urlopen(_url) as _f:
            _raw = json.load(_f)
    else:
        with open(_url) as _f:
            _raw = json.load(_f)
    # Extract accessions before converting site keys
    nt_accessions = _raw.pop("_accessions", {})
    # Convert string site keys to int
    reference_nt_table = {
        key: {int(k): v for k, v in sites.items()}
        for key, sites in _raw.items()
    }
    return nt_accessions, reference_nt_table


@app.cell
def _(nt_accessions, reference_nt_table, segment, subtype):
    _ref_subtype = subtype.value if segment.value in ("HA", "NA") else "all"
    _key = f"{segment.value}:{_ref_subtype}"
    reference_nt = reference_nt_table.get(_key, {})
    reference_nt_accession = nt_accessions.get(_key, "")
    return reference_nt, reference_nt_accession


@app.cell
def _(df, host, min_count, mut_class, reference_nt, segment, subtype):
    _df = df[df["segment"] == segment.value].copy()

    # Host filter
    _df = _df[_df["host"] == host.value]

    # Subtype filter
    if segment.value in ("HA", "NA"):
        _df = _df[_df["subtype"] == subtype.value]
    else:
        _df = _df[_df["subtype"] == "all"]

    # Mutation class filter
    if mut_class.value != "all":
        _df = _df[_df["mut_class"] == mut_class.value]

    # Count filter
    _count_mask = (_df["expected_count"] >= min_count.value) | (_df["actual_count"] >= min_count.value)
    _df = _df[_count_mask]

    # Reference filter: only show mutations where wt_nt matches the reference sequence
    if reference_nt:
        _df = _df[_df["site"].map(reference_nt) == _df["wt_nt"]]

    plot_data = _df.copy()
    return (plot_data,)


@app.cell
def _(alt, plot_data):
    NT_ORDER = ["A", "C", "G", "T"]

    color_scale = alt.Scale(
        scheme="redblue",
        domainMid=0,
        domainMin=-6,
        domainMax=2,
        clamp=True,
    )

    # Brush on quantitative x so numeric interval selection works correctly.
    # Set an initial value so the heatmap doesn't try to render all sites at once —
    # with site:O and width=Step(16), a 2000+ site segment would be ~36 000 px wide.
    _all_sites = sorted(plot_data["site"].unique()) if len(plot_data) > 0 else [1, 100]
    _init_brush_max = _all_sites[min(199, len(_all_sites) - 1)]
    site_brush = alt.selection_interval(
        name="site_brush",
        encodings=["x"],
        value={"x": [_all_sites[0], _init_brush_max]},
        mark=alt.BrushConfig(stroke="black", strokeWidth=2),
    )

    # --- Panel 1: site zoom bar ---
    _sites_df = plot_data[["site"]].drop_duplicates().sort_values("site")
    site_zoom_bar = (
        alt.Chart(_sites_df)
        .mark_rect(color="#BB86FC", opacity=0.6)
        .encode(x=alt.X("site:Q", title=None, axis=None, scale=alt.Scale(zero=False, nice=False)))
        .add_params(site_brush)
        .properties(width=900, height=20, title="Drag to zoom")
    )

    # --- Panel 2: site-level line plot (mean fitness per site) ---
    _site_agg = (
        plot_data
        .groupby("site", as_index=False)
        .agg(mean_fitness=("delta_fitness", "mean"), gene=("gene", lambda x: "/".join(sorted(x.dropna().unique()))))
    )
    lineplot = (
        alt.Chart(_site_agg)
        .mark_line(point=alt.OverlayMarkDef(color="#BB86FC", filled=True), opacity=0.9, color="#BB86FC")
        .transform_filter(site_brush)
        .encode(
            x=alt.X("site:Q", title="Site", scale=alt.Scale(zero=False, nice=False)),
            y=alt.Y("mean_fitness:Q", title="Mean fitness effect",
                    scale=alt.Scale(zero=False)),
            tooltip=[
                alt.Tooltip("site:Q", title="Site"),
                alt.Tooltip("gene:N", title="Gene"),
                alt.Tooltip("mean_fitness:Q", title="Mean fitness effect", format=".3f"),
            ],
        )
        .properties(width=900, height=130)
    )

    # --- Panel 3: per-mutation heatmap ---
    # Background: grey bar per site to indicate sites that have data
    _sites_only = plot_data[["site"]].drop_duplicates()
    background = (
        alt.Chart(_sites_only)
        .mark_bar(color="lightgray", opacity=0.2, width={"band": 1.0})
        .transform_filter(site_brush)
        .encode(x=alt.X("site:O", title="Site",
                         axis=alt.Axis(labelOverlap=False, labelAngle=-90)))
    )

    # Colored rectangles
    heatmap_rects = (
        alt.Chart(plot_data)
        .mark_rect()
        .transform_filter(site_brush)
        .encode(
            x=alt.X("site:O", title="Site",
                    axis=alt.Axis(labelOverlap=False, labelAngle=-90)),
            y=alt.Y("mut_nt:N", title="Mutant NT",
                    scale=alt.Scale(domain=NT_ORDER),
                    axis=alt.Axis(grid=True, tickBand="extent",
                                  gridColor="#333333", gridWidth=2)),
            color=alt.Color("delta_fitness:Q", scale=color_scale, title="Fitness effect"),
            tooltip=[
                alt.Tooltip("site:Q", title="Site"),
                alt.Tooltip("gene:N", title="Gene"),
                alt.Tooltip("wt_nt:N", title="Ref nt"),
                alt.Tooltip("mut_nt:N", title="Mut nt"),
                alt.Tooltip("mut_class:N", title="Mutation class"),
                alt.Tooltip("delta_fitness:Q", title="Fitness effect", format=".3f"),
                alt.Tooltip("actual_count:Q", title="Actual count", format=".1f"),
                alt.Tooltip("expected_count:Q", title="Expected count", format=".1f"),
            ],
        )
    )

    # Wildtype (reference) markers: "x" at each (site, wt_nt) position
    _wt_df = plot_data[["site", "wt_nt"]].drop_duplicates()
    wt_marks = (
        alt.Chart(_wt_df)
        .mark_text(text="x", color="black", fontSize=10)
        .transform_filter(site_brush)
        .encode(
            x=alt.X("site:O", axis=alt.Axis(labelOverlap=False, labelAngle=-90)),
            y=alt.Y("wt_nt:N", scale=alt.Scale(domain=NT_ORDER)),
        )
    )

    heatmap = (background + heatmap_rects + wt_marks).properties(width=alt.Step(16), height=70)

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
    host,
    min_count,
    mo,
    mut_class,
    plot_data,
    reference_nt_accession,
    segment,
    subtype,
):
    _n_muts = len(plot_data)
    _n_sites = plot_data["site"].nunique()
    _acc = reference_nt_accession or "—"
    summary = mo.md(f"**{_n_muts} mutations shown across {_n_sites} sites (reference sequence: {_acc}).**")

    if segment.value in ("HA", "NA"):
        controls = mo.hstack(
            [segment, subtype, host, mut_class, min_count],
            gap=2,
        )
    else:
        controls = mo.hstack(
            [segment, host, mut_class, min_count],
            gap=2,
        )

    description = mo.md("""
    ## Fitness effects of nucleotide-level mutations to influenza virus protein-coding genes

    This dashboard shows the fitness effects of single-nucleotide mutations to influenza
    protein-coding sequences, as estimated in Haddox et al., 2026.

    **Fitness effects** are estimated as the log ratio of actual to expected mutation counts
    at each site. Actual counts correspond to the number of times a mutation is observed
    to occur along the branches of a phylogenetic tree, and expected counts are derived from a
    model that estimates how many times a mutation is expected to occur under neutrality.
    Negative values (blue) indicate the mutation is deleterious, while positive values (red)
    indicate the mutation is beneficial.

    **How to use:**
    - Select a **segment**, **subtype**, and **host** with the dropdowns.
    - Use **Mutation class** to focus on synonymous (default), nonsynonymous, or all mutations.
    - Use **Min count** to hide sites with low counts. The dashboard will only show mutations
      with at least the indicated number of actual or expected counts.
    - **Drag** the purple zoom bar to pan and zoom into a region of interest.
    - The **line plot** shows the mean fitness effect of all mutations at a site.
    - The **heatmap** shows site-specific fitness effects of mutations. An **×** marks
      the reference nucleotide at each site.
    - **Hover** over any point or cell for detailed values.

    The heatmap shows mutations to a specific reference sequence, the accession of which is
    given in the summary line above the plot. The full dataset includes additional mutations.
    """)

    mo.vstack([controls, summary, mo.ui.altair_chart(chart), description])
    return


if __name__ == "__main__":
    app.run()
