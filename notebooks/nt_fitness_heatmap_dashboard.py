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
    HOSTS = ["all", "human", "avian"]
    host = mo.ui.dropdown(options=HOSTS, value="all", label="Host")
    MUT_CLASSES = ["all", "nonsynonymous", "synonymous"]
    mut_class = mo.ui.dropdown(options=MUT_CLASSES, value="all", label="Mutation class")
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
    # Convert string site keys to int
    reference_nt_table = {
        key: {int(k): v for k, v in sites.items()}
        for key, sites in _raw.items()
    }
    return (reference_nt_table,)


@app.cell
def _(protein, reference_nt_table, subtype):
    _GENE_TO_SEGMENT = {
        "HA": "HA", "NA": "NA", "NP": "NP",
        "PA": "PA", "PB1": "PB1", "PB2": "PB2",
        "M1": "MP", "M2": "MP",
        "NS1": "NS", "NEP": "NS",
    }
    _ref_subtype = subtype.value if protein.value in ("HA", "NA") else "all"
    _segment = _GENE_TO_SEGMENT[protein.value]
    reference_nt = reference_nt_table.get(f"{_segment}:{_ref_subtype}", {})
    return (reference_nt,)


@app.cell
def _(df, host, min_count, mut_class, protein, reference_nt, subtype):
    _df = df[df["gene"] == protein.value].copy()

    # Host filter
    _df = _df[_df["host"] == host.value]

    # Subtype filter
    if protein.value in ("HA", "NA"):
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
        .groupby("site", as_index=False)["delta_fitness"]
        .mean()
        .rename(columns={"delta_fitness": "mean_fitness"})
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
def _(chart, host, min_count, mo, mut_class, plot_data, protein, subtype):
    _n_muts = len(plot_data)
    _n_sites = plot_data["site"].nunique()
    summary = mo.md(f"**{_n_muts} mutations shown across {_n_sites} sites.**")

    if protein.value in ("HA", "NA"):
        controls = mo.hstack(
            [protein, subtype, host, mut_class, min_count],
            gap=2,
        )
    else:
        controls = mo.hstack(
            [protein, host, mut_class, min_count],
            gap=2,
        )

    mo.vstack([controls, summary, mo.ui.altair_chart(chart)])
    return


if __name__ == "__main__":
    app.run()
