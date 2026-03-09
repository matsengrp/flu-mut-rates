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
    # Explode overlapping-ORF rows (e.g. "NEP;NS1") on all four joined columns
    # Read via urllib + StringIO to prevent pandas double-decompressing gzip responses
    _url = f"{DATA_DIR}/aa_fitness_effects.csv"
    if _url.startswith("http"):
        with _urllib.urlopen(_url) as _r:
            df = pd.read_csv(io.StringIO(_r.read().decode("utf-8")), keep_default_na=False)
    else:
        df = pd.read_csv(_url, keep_default_na=False)
    for _col in ["gene", "codon_site", "wt_aa", "mut_aa"]:
        df[_col] = df[_col].str.split(";")
    df = df.explode(["gene", "codon_site", "wt_aa", "mut_aa"])
    df["codon_site"] = df["codon_site"].astype(int)
    return (df,)


@app.cell
def _(DATA_DIR):
    import json
    import urllib.request
    _url = f"{DATA_DIR}/reference_aa.json"
    if _url.startswith("http"):
        with urllib.request.urlopen(_url) as _f:
            _raw = json.load(_f)
    else:
        with open(_url) as _f:
            # Keys are "GENE:SUBTYPE"; codon sites are string keys -> convert to int
            _raw = json.load(_f)
    reference_aa_table = {
        key: {int(k): v for k, v in sites.items()}
        for key, sites in _raw.items()
    }
    return (reference_aa_table,)


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
    min_count = mo.ui.slider(
        start=0, stop=200, step=5, value=10,
        label="Min count (actual or expected)", show_value=True,
    )
    return host, min_count


@app.cell
def _(protein, reference_aa_table, subtype):
    _ref_subtype = subtype.value if protein.value in ("HA", "NA") else "all"
    reference_aa = reference_aa_table.get(f"{protein.value}:{_ref_subtype}", {})
    return (reference_aa,)


@app.cell
def _(df, host, min_count, protein, reference_aa, subtype):
    _df = df[df["gene"] == protein.value].copy()

    # Host filter
    _df = _df[_df["host"] == host.value]

    # Subtype filter
    if protein.value in ("HA", "NA"):
        _df = _df[_df["subtype"] == subtype.value]
    else:
        _df = _df[_df["subtype"] == "all"]

    # Count filter
    _count_mask = (_df["expected_count"] >= min_count.value) | (_df["actual_count"] >= min_count.value)
    _df = _df[_count_mask]

    # Reference filter: only show mutations where wt_aa matches the reference
    # sequence amino acid (applied to all mutation classes)
    if reference_aa:
        _df = _df[_df["codon_site"].map(reference_aa) == _df["wt_aa"]]

    plot_data = _df.copy()
    return (plot_data,)


@app.cell
def _(alt, plot_data):
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
        .encode(x=alt.X("codon_site:Q", title=None, axis=None, scale=alt.Scale(zero=False, nice=False)))
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
            x=alt.X("codon_site:Q", title="Site", scale=alt.Scale(zero=False, nice=False)),
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

    # Colored rectangles — color scale clamps display; actual values shown in tooltip
    heatmap_rects = (
        alt.Chart(plot_data)
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
def _(chart, host, min_count, mo, plot_data, protein, subtype):
    _n_muts = len(plot_data)
    _n_sites = plot_data["codon_site"].nunique()
    summary = mo.md(f"**{_n_muts} mutations shown across {_n_sites} sites.**")

    if protein.value in ("HA", "NA"):
        controls = mo.hstack(
            [protein, subtype, host, min_count],
            gap=2,
        )
    else:
        controls = mo.hstack(
            [protein, host, min_count],
            gap=2,
        )

    mo.vstack([controls, summary, mo.ui.altair_chart(chart)])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
