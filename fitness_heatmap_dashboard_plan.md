# Plan: Interactive Fitness Effects Heatmap Dashboard (Marimo Notebook)

## Goal

Build a marimo notebook (`notebooks/fitness_heatmap_dashboard.py`) that provides an
interactive dashboard for exploring amino-acid fitness effects across influenza proteins.
The dashboard will show a per-site, per-mutation heatmap of `delta_fitness` values,
with UI controls to select a protein and filter by count thresholds.

---

## Data Overview

**Source:** `results/aa_fitness_effects.csv`

**Columns:**
- `host` — `all`, `human`, or `avian`
- `subtype` — e.g. `H1`, `H3`, `N1`, `all`
- `segment` — e.g. `HA`, `NA`, `PB2`, ...
- `gene` — e.g. `HA`, `NA`, `NP`, `M1`, `M2`, `NS1`, `NEP`, `PA`, `PB1`, `PB2`
  - Some rows have multi-gene entries like `M2;M1` or `NEP;NS1` for mutations in overlapping ORFs
- `codon_site` — amino-acid position (can also be `;`-joined for overlapping ORFs)
- `wt_aa` — wildtype amino acid (can also be `;`-joined)
- `mut_aa` — mutant amino acid (can also be `;`-joined)
- `aa_mut` — string like `M1I`
- `mut_class` — `synonymous`, `nonsynonymous`, or `nonsense`
- `actual_count` — observed mutation count
- `expected_count` — expected count from neutral model
- `delta_fitness` — log fitness effect

**Key data issues to handle:**
- Overlapping ORF rows (e.g. `gene = "NEP;NS1"`) must be exploded on
  `gene`, `codon_site`, `wt_aa`, `mut_aa` simultaneously, then filtered
  to the gene of interest.
- Only show mutations where `mut_aa != wt_aa` (i.e. exclude synonymous identity rows;
  these arise from the overlapping ORF explosion).
- For HA and NA, there are multiple subtypes; the user may want to select a subtype.
  Default to `H1` for HA and `N1` for NA; internal segments use `all`.
- For the heatmap, focus on `host == "all"` and `mut_class == "nonsynonymous"`.

---

## Dashboard Layout

```
+--------------------------------------------------+
| [Protein selector dropdown]                      |
| [Subtype selector] (only shown for HA/NA)        |
| [Min expected count slider]                      |
| [Min actual count slider / "either" toggle]      |
+--------------------------------------------------+
| Site zoom bar (brush to select site range)       |
+--------------------------------------------------+
| Line plot: site (x) vs. mean delta_fitness (y)  |
+--------------------------------------------------+
| Heatmap: site (x-axis) x mutant AA (y-axis)     |
|          colored by delta_fitness                |
+--------------------------------------------------+
| Summary stats (N mutations shown, site range)    |
+--------------------------------------------------+
```

All three chart panels (zoom bar, line plot, heatmap) share the same x-axis and
are composed with `alt.vconcat`. Brushing on the zoom bar zooms the x-axis of the
line plot and heatmap.

---

## Step-by-Step Implementation Plan

### Step 1: Data Loading and Preprocessing

1. Load `results/aa_fitness_effects.csv` with `keep_default_na=False`.
2. Filter to `host == "all"` and `mut_class == "nonsynonymous"`.
3. Handle overlapping ORF rows:
   - Detect rows where `gene` contains `;`.
   - For those rows, split and explode `gene`, `codon_site`, `wt_aa`, `mut_aa`
     simultaneously (as done in `analyze_fitness_effects.ipynb`).
   - After exploding, remove rows where `wt_aa == mut_aa` (synonymous in the
     exploded gene — these are nonsynonymous in the other overlapping gene).
4. Cast `codon_site` to `int`.
5. Normalize `gene` values to uppercase (raw data has mixed case, e.g. `ha` vs `HA`).
6. Remove entries for `PB1;PB1-F2` (as done in existing notebook).
7. Store the preprocessed full DataFrame as a module-level variable (computed once).

**Available proteins after preprocessing:**
`HA`, `NA`, `NP`, `PA`, `PB1`, `PB2`, `M1`, `M2`, `NS1`, `NEP`

### Step 2: UI Controls (Marimo Reactive Widgets)

```python
import marimo as mo

protein = mo.ui.dropdown(
    options=["HA", "NA", "NP", "PA", "PB1", "PB2", "M1", "M2", "NS1", "NEP"],
    value="HA",
    label="Protein",
)

# Subtype dropdown — only meaningful for HA and NA
# HA subtypes: H1, H3, H5, H7, H9
# NA subtypes: N1, N2, N6, N8, N9
subtype = mo.ui.dropdown(
    options=["H1"],   # dynamically updated based on protein
    value="H1",
    label="Subtype",
)

min_expected_count = mo.ui.slider(
    start=0, stop=200, step=5, value=10,
    label="Min expected count",
)

min_actual_count = mo.ui.slider(
    start=0, stop=200, step=5, value=0,
    label="Min actual count",
)

count_filter_mode = mo.ui.radio(
    options=["either (expected OR actual)", "both (expected AND actual)"],
    value="either (expected OR actual)",
    label="Count filter mode",
)
```

- The subtype dropdown should only be shown/active when the protein is `HA` or `NA`.
- When protein changes to/from HA/NA, the subtype dropdown options update reactively.

### Step 3: Data Subsetting (Reactive)

Given the user's selections, subset the preprocessed DataFrame:

```python
def get_plot_data(df, protein, subtype, min_expected, min_actual, filter_mode):
    data = df[df["gene"] == protein].copy()

    # Apply subtype filter for HA/NA
    if protein in ("HA", "NA"):
        data = data[data["subtype"] == subtype]
    else:
        data = data[data["subtype"] == "all"]

    # Apply count filter
    if filter_mode == "either":
        mask = (data["expected_count"] >= min_expected) | (data["actual_count"] >= min_actual)
    else:
        mask = (data["expected_count"] >= min_expected) & (data["actual_count"] >= min_actual)
    data = data[mask]

    return data
```

### Step 4: Combined Line Plot + Heatmap with Altair

The visualization is a three-panel `alt.vconcat` chart, directly modeled on the
`_lineplot_and_heatmap` function from
[multidms/plot.py](https://github.com/matsengrp/multidms/blob/cd2d73214f8163cc9c00e1ee30383bc7f65fed22/multidms/plot.py#L49).

#### Shared setup

```python
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY*")

# Diverging color scale clipped to [-6, 2], centered at 0
color_scale = alt.Scale(
    scheme="redblue",
    domainMid=0,
    domainMin=-6,
    domainMax=2,
    clamp=True,
)

# Brush for site-range zoom (defined once, shared across all panels)
site_brush = alt.selection_interval(
    encodings=["x"],
    mark=alt.BrushConfig(stroke="black", strokeWidth=2),
)

# Shared x-encoding that filters to the brushed range
x_zoom = alt.X(
    "codon_site:O",
    scale=alt.Scale(domain={"selection": site_brush.name, "encoding": "x"}),
    title="Site",
    axis=alt.Axis(labelOverlap=True),
)
```

#### Panel 1: Site zoom bar

A thin bar chart (one rectangle per site) that the user brushes to zoom into a
site range. All sites are shown regardless of count filters so the user has
full context of the protein length.

```python
site_zoom_bar = (
    alt.Chart(data[["codon_site"]].drop_duplicates())
    .mark_rect(color="steelblue", opacity=0.6)
    .encode(
        x=alt.X("codon_site:O", title=None, axis=None),
    )
    .add_params(site_brush)
    .properties(width=800, height=20, title="drag to zoom")
)
```

#### Panel 2: Site-level line plot

Aggregates per-mutation `delta_fitness` values to a per-site mean (or median),
then plots a point + line chart. X-axis is linked to the zoom brush.

```python
lineplot = (
    alt.Chart(data)
    .mark_line(point=True, opacity=0.7)
    .encode(
        x=x_zoom,
        y=alt.Y(
            "mean(delta_fitness):Q",
            title="Mean fitness effect",
            scale=alt.Scale(zero=False),
        ),
        tooltip=[
            alt.Tooltip("codon_site:O", title="Site"),
            alt.Tooltip("mean(delta_fitness):Q", title="Mean fitness effect", format=".3f"),
            alt.Tooltip("count():Q", title="N mutations"),
        ],
    )
    .properties(width=800, height=120)
)
```

#### Panel 3: Heatmap

Each cell is a (site, mutant AA) pair colored by `delta_fitness`. Three layers
are composed with `+`:

1. **Background layer** — gray rectangles for all (site, AA) combinations that
   *could* exist (cross-join of sites × AA_ORDER), indicating missing/filtered data.
2. **Data layer** — colored rectangles for mutations that pass the count filter.
3. **Wildtype markers** — black `x` text at (site, wt_aa) positions.

```python
# Background: all possible (site, AA) positions
all_positions = pd.MultiIndex.from_product(
    [data["codon_site"].unique(), AA_ORDER], names=["codon_site", "mut_aa"]
).to_frame(index=False)

background = (
    alt.Chart(all_positions)
    .mark_rect(color="lightgray", opacity=0.3)
    .encode(
        x=x_zoom,
        y=alt.Y("mut_aa:N", sort=AA_ORDER, title="Mutant AA"),
    )
)

# Data layer
heatmap_rects = (
    alt.Chart(data)
    .mark_rect()
    .encode(
        x=x_zoom,
        y=alt.Y("mut_aa:N", sort=AA_ORDER, title="Mutant AA"),
        color=alt.Color("delta_fitness:Q", scale=color_scale, title="Fitness effect"),
        tooltip=[
            alt.Tooltip("codon_site:O", title="Site"),
            alt.Tooltip("wt_aa:N", title="WT AA"),
            alt.Tooltip("mut_aa:N", title="Mut AA"),
            alt.Tooltip("delta_fitness:Q", title="Delta fitness", format=".3f"),
            alt.Tooltip("actual_count:Q", title="Actual count", format=".1f"),
            alt.Tooltip("expected_count:Q", title="Expected count", format=".1f"),
        ],
    )
)

# Wildtype markers
wt_df = data[["codon_site", "wt_aa"]].drop_duplicates()
wt_marks = (
    alt.Chart(wt_df)
    .mark_text(text="x", color="black", fontSize=8)
    .encode(
        x=x_zoom,
        y=alt.Y("wt_aa:N", sort=AA_ORDER),
    )
)

heatmap = (background + heatmap_rects + wt_marks).properties(width=800, height=300)
```

#### Final composition

```python
chart = (
    alt.vconcat(site_zoom_bar, lineplot, heatmap, spacing=5)
    .resolve_scale(color="independent", x="shared")
    .configure_axis(grid=False, labelOverlap="parity")
    .configure(padding=10)
)
```

The `x="shared"` resolution ensures the zoom brush on the zoom bar applies to
both the line plot and the heatmap simultaneously.

**Handling large proteins:** HA has ~560 sites. With `width=800` pixels in a
fixed-width layout the sites will be compressed; the zoom bar lets users drill
into specific regions. Use `alt.Step(12)` for heatmap cell width as an alternative
if per-site readability at full zoom is a priority.

### Step 5: Dashboard Assembly in Marimo

```python
@mo.cell
def dashboard():
    plot_data = get_plot_data(
        preprocessed_df,
        protein.value,
        subtype.value,
        min_expected_count.value,
        min_actual_count.value,
        count_filter_mode.value,
    )
    chart = make_lineplot_and_heatmap(plot_data)

    n_muts = len(plot_data)
    n_sites = plot_data["codon_site"].nunique()
    summary = mo.md(f"**Showing {n_muts} mutations across {n_sites} sites.**")

    # Show subtype selector only for HA/NA
    if protein.value in ("HA", "NA"):
        controls = mo.vstack([protein, subtype, min_expected_count, min_actual_count, count_filter_mode])
    else:
        controls = mo.vstack([protein, min_expected_count, min_actual_count, count_filter_mode])

    return mo.vstack([controls, mo.altair_chart(chart), summary])
```

### Step 6: Subtype Selector Reactivity

Because marimo reactive cells re-run when dependencies change, the subtype dropdown
options need to update when the protein changes. This can be done by defining a
reactive cell that recomputes the subtype dropdown:

```python
@mo.cell
def subtype_options(protein):
    HA_SUBTYPES = ["H1", "H3", "H5", "H7", "H9"]
    NA_SUBTYPES = ["N1", "N2", "N6", "N8", "N9"]
    if protein.value == "HA":
        opts = HA_SUBTYPES
        default = "H1"
    elif protein.value == "NA":
        opts = NA_SUBTYPES
        default = "N1"
    else:
        opts = ["all"]
        default = "all"
    return mo.ui.dropdown(options=opts, value=default, label="Subtype")
```

---

## File Structure

```
notebooks/
└── fitness_heatmap_dashboard.py   # New marimo notebook
```

Run with:
```bash
marimo edit notebooks/fitness_heatmap_dashboard.py
# or for a read-only app:
marimo run notebooks/fitness_heatmap_dashboard.py
```

---

## Dependencies

All dependencies are already in the `flu-syn-rates` conda environment:
- `marimo` (check; add to `environment.yml` if missing)
- `altair`
- `pandas`
- `numpy`

---

## Open Questions / Design Decisions

1. **Site-range filtering:** For long proteins (HA, PB2, PB1), the heatmap will be
   very wide. Should we add a site-range slider (start/end site) to zoom into a region?
   Recommended: yes, add as a secondary control.

2. **Color scale clipping:** The current notebook clips `delta_fitness` to `[-7, 2]`.
   Should the heatmap use the same clip, or show a continuous scale with saturation?
   Recommend: clip to `[-6, 2]` and use a diverging scale.

3. **Grey-out vs. hide:** Mutations that don't pass the count threshold — should they
   be shown as grey (indicating "no data") or hidden entirely?
   Recommend: hide (simpler); grey-out is nicer but requires a background layer.

4. **NA segment:** Some NA subtypes are present in the data. The `gene` column stores
   `NA` (uppercase) and the segment column is also `NA`. Need to handle the `na` (lowercase)
   issue (raw data has mixed case) via normalization in Step 1.

5. **Subtype selection for internal segments:** These all use `subtype == "all"` so no
   subtype UI is needed.

6. **`host` selection:** Currently plan to fix `host == "all"`. Could add a `host`
   dropdown (`all`, `human`, `avian`) as a future enhancement.
