"""Render single-panel fitness-effect heatmap PNGs for the manuscript figures.

Mirrors the heatmap panel of notebooks/{aa,nt}_fitness_heatmap_dashboard.py but
trimmed: no zoom bar, no line-plot panel, fixed site range applied as a pandas
filter. Outputs feed notebooks/compose_figures.ipynb.

Visual specs (color scale, AA/NT order, tick styling) are duplicated from the
dashboards. Keep in sync with notebooks/{aa,nt}_fitness_heatmap_dashboard.py.
"""

import argparse
import json
from pathlib import Path

import altair as alt
import pandas as pd
import vl_convert as vlc


AA_ORDER = [
    "K", "R", "H",
    "D", "E",
    "N", "Q", "S", "T",
    "P", "G", "A", "V", "L", "I", "M",
    "F", "W", "Y",
    "C",
    "*",
]

NT_ORDER = ["A", "C", "G", "T"]

COLOR_SCALE = alt.Scale(
    scheme="redblue",
    domainMid=0,
    domainMin=-6,
    domainMax=2,
    clamp=True,
)

# (site_min, site_max) inclusive. nt_heatmap defaults to the dashboard's default
# view (HA:H1, synonymous); update if a different segment is wanted.
SPECS = {
    "aa_heatmap": {
        "level": "aa", "gene": "PA", "subtype": "all",
        "host": "all", "min_count": 10, "sites": (506, 540),
    },
    "nt_heatmap": {
        "level": "nt", "segment": "MP", "subtype": "all",
        "host": "all", "mut_class": "synonymous",
        "min_count": 10, "sites": (51, 213),
    },
    "MP_5SS": {
        "level": "nt", "segment": "MP", "subtype": "all",
        "host": "all", "mut_class": "synonymous",
        "min_count": 10, "sites": (10, 42),
    },
    "MP_3SS": {
        "level": "nt", "segment": "MP", "subtype": "all",
        "host": "all", "mut_class": "synonymous",
        "min_count": 10, "sites": (696, 763),
    },
    "NS_5SS": {
        "level": "nt", "segment": "NS", "subtype": "all",
        "host": "all", "mut_class": "synonymous",
        "min_count": 10, "sites": (9, 60),
    },
    "NS_3SS": {
        "level": "nt", "segment": "NS", "subtype": "all",
        "host": "all", "mut_class": "synonymous",
        "min_count": 10, "sites": (468, 501),
    },
}


def _load_reference(path: Path) -> dict:
    with open(path) as f:
        raw = json.load(f)
    raw.pop("_accessions", None)
    return {key: {int(k): v for k, v in sites.items()} for key, sites in raw.items()}


def _load_aa_df(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, keep_default_na=False)
    # Explode overlapping-ORF rows (e.g. "NEP;NS1") on the four joined columns
    for col in ["gene", "codon_site", "wt_aa", "mut_aa"]:
        df[col] = df[col].str.split(";")
    df = df.explode(["gene", "codon_site", "wt_aa", "mut_aa"])
    df["codon_site"] = df["codon_site"].astype(int)
    return df


def _filter_aa(df: pd.DataFrame, spec: dict, reference_aa: dict) -> pd.DataFrame:
    out = df[df["gene"] == spec["gene"]]
    out = out[out["host"] == spec["host"]]
    sub_filter = spec["subtype"] if spec["gene"] in ("HA", "NA") else "all"
    out = out[out["subtype"] == sub_filter]
    mc = spec["min_count"]
    out = out[(out["expected_count"] >= mc) | (out["actual_count"] >= mc)]
    if reference_aa:
        out = out[out["codon_site"].map(reference_aa) == out["wt_aa"]]
    lo, hi = spec["sites"]
    out = out[(out["codon_site"] >= lo) & (out["codon_site"] <= hi)]
    return out


def _filter_nt(df: pd.DataFrame, spec: dict, reference_nt: dict) -> pd.DataFrame:
    out = df[df["segment"] == spec["segment"]]
    out = out[out["host"] == spec["host"]]
    sub_filter = spec["subtype"] if spec["segment"] in ("HA", "NA") else "all"
    out = out[out["subtype"] == sub_filter]
    if spec["mut_class"] != "all":
        out = out[out["mut_class"] == spec["mut_class"]]
    mc = spec["min_count"]
    out = out[(out["expected_count"] >= mc) | (out["actual_count"] >= mc)]
    if reference_nt:
        out = out[out["site"].map(reference_nt) == out["wt_nt"]]
    lo, hi = spec["sites"]
    out = out[(out["site"] >= lo) & (out["site"] <= hi)]
    return out


def _build_chart(plot_data: pd.DataFrame, *, level: str, height: int) -> alt.Chart:
    if level == "aa":
        x_field, y_field, wt_field = "codon_site", "mut_aa", "wt_aa"
        y_domain, y_title = AA_ORDER, "Mutant AA"
        label_angle = 0
        label_overlap = True
    else:
        x_field, y_field, wt_field = "site", "mut_nt", "wt_nt"
        y_domain, y_title = NT_ORDER, "Mutant NT"
        label_angle = -90
        label_overlap = False

    x_axis = alt.Axis(labelOverlap=label_overlap, labelAngle=label_angle)

    sites_only = plot_data[[x_field]].drop_duplicates()
    background = (
        alt.Chart(sites_only)
        .mark_bar(color="lightgray", opacity=0.2, width={"band": 1.0})
        .encode(x=alt.X(f"{x_field}:O", title="Site", axis=x_axis))
    )

    heatmap_rects = (
        alt.Chart(plot_data)
        .mark_rect()
        .encode(
            x=alt.X(f"{x_field}:O", title="Site", axis=x_axis),
            y=alt.Y(
                f"{y_field}:N",
                title=y_title,
                scale=alt.Scale(domain=y_domain),
                axis=alt.Axis(grid=True, tickBand="extent",
                              gridColor="#333333", gridWidth=2),
            ),
            color=alt.Color("delta_fitness:Q", scale=COLOR_SCALE,
                            title="Fitness effect"),
        )
    )

    wt_df = plot_data[[x_field, wt_field]].drop_duplicates()
    wt_marks = (
        alt.Chart(wt_df)
        .mark_text(text="x", color="black", fontSize=10)
        .encode(
            x=alt.X(f"{x_field}:O", axis=x_axis),
            y=alt.Y(f"{wt_field}:N", scale=alt.Scale(domain=y_domain)),
        )
    )

    # Order matters: a bare .configure(...) call replaces the entire config
    # object (wiping out configure_axis / configure_legend), so put all the
    # subconfigs together via configure(...) once.
    return (
        (background + heatmap_rects + wt_marks)
        .properties(width=alt.Step(16), height=height)
        .configure(
            background="#2b2b2b",
            padding=10,
            view=alt.ViewConfig(stroke=None),
            axis=alt.AxisConfig(
                grid=False, labelOverlap="parity",
                domainColor="white", tickColor="white",
                labelColor="white", titleColor="white",
                labelFontWeight="bold", titleFontWeight="bold",
            ),
            legend=alt.LegendConfig(
                labelColor="white", titleColor="white",
                labelFontWeight="bold", titleFontWeight="bold",
            ),
        )
    )


def render(name: str, *, aa_df, nt_df, ref_aa, ref_nt, out_dir: Path) -> Path:
    spec = SPECS[name]
    if spec["level"] == "aa":
        sub = spec["subtype"] if spec["gene"] in ("HA", "NA") else "all"
        ref = ref_aa.get(f"{spec['gene']}:{sub}", {})
        data = _filter_aa(aa_df, spec, ref)
        height = 320
    else:
        sub = spec["subtype"] if spec["segment"] in ("HA", "NA") else "all"
        ref = ref_nt.get(f"{spec['segment']}:{sub}", {})
        data = _filter_nt(nt_df, spec, ref)
        height = 70

    if len(data) == 0:
        raise SystemExit(
            f"{name}: no rows after filtering — check spec {spec}"
        )

    chart = _build_chart(data, level=spec["level"], height=height)
    png_bytes = vlc.vegalite_to_png(chart.to_json(), scale=2)
    out_path = out_dir / f"{name}.png"
    out_path.write_bytes(png_bytes)
    print(f"  -> {out_path}  ({len(data)} rows, {data[data.columns[4]].nunique()} sites)")
    return out_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--aa-csv", default="results/aa_fitness_effects.csv")
    p.add_argument("--nt-csv", default="results/nt_fitness_effects.csv")
    p.add_argument("--aa-ref", default="results/reference_aa.json")
    p.add_argument("--nt-ref", default="results/reference_nt.json")
    p.add_argument("--output-dir", default="results/figures")
    p.add_argument("--name", action="append", default=None,
                   help="Render only the named snapshot(s). Repeatable. "
                        "Default: render all.")
    args = p.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    alt.data_transformers.disable_max_rows()

    names = args.name or list(SPECS.keys())
    unknown = [n for n in names if n not in SPECS]
    if unknown:
        raise SystemExit(f"Unknown snapshot name(s): {unknown}. "
                         f"Available: {list(SPECS)}")

    # Lazy-load CSVs only if needed
    need_aa = any(SPECS[n]["level"] == "aa" for n in names)
    need_nt = any(SPECS[n]["level"] == "nt" for n in names)
    aa_df = _load_aa_df(Path(args.aa_csv)) if need_aa else None
    nt_df = pd.read_csv(args.nt_csv, keep_default_na=False) if need_nt else None
    ref_aa = _load_reference(Path(args.aa_ref)) if need_aa else {}
    ref_nt = _load_reference(Path(args.nt_ref)) if need_nt else {}

    for name in names:
        print(f"Rendering {name}...")
        render(name, aa_df=aa_df, nt_df=nt_df, ref_aa=ref_aa, ref_nt=ref_nt,
               out_dir=out_dir)


if __name__ == "__main__":
    main()
