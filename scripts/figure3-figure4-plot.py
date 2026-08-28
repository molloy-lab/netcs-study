#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.2,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

METHOD_FAMILY_COLORS = {
    "ASTRAL": "#666666",
    "CAMUS": "#4477AA",
    "NANQUPLUS": "#CC6677",
    "NET-CS": "#117733",
}

METHOD_COLORS = {
    "ASTRAL": "#777777",
    "CAMUS": "#4477AA",

    # Estimated ToBs: lighter
    "NANQUPLUS-estimated-tob": "#DC8FA0",
    "NET-CS-estimated-tob": "#78AD87",

    # True ToBs: darker
    "NANQUPLUS-true-tob": "#8F2946",
    "NET-CS-true-tob": "#075A2A",
}

FAMILY_LEGEND_LABELS = {
    "ASTRAL": "ASTRAL",
    "CAMUS": "CAMUS",
    "NANQUPLUS": "NQ+",
    "NET-CS": "NetCS",
}

SOURCE_LEGEND_COLORS = {
    "estimated-tob": "#B5B5B5",
    "true-tob": "#3F3F3F",
}

SOURCE_LEGEND_LABELS = {
    "estimated-tob": "Estimated ToBs",
    "true-tob": "True ToBs",
}

MISSING_GROUP_MARKER = "NR"

GROUP_GAP = 1.45

HWCD_METHOD_ORDER = [
    "ASTRAL",
    "CAMUS",
    "NANQUPLUS-estimated-tob",
    "NET-CS-estimated-tob",
    "NANQUPLUS-true-tob",
    "NET-CS-true-tob",
]

HYBRID_METHOD_ORDER = [
    "CAMUS",
    "NANQUPLUS-estimated-tob",
    "NET-CS-estimated-tob",
    "NANQUPLUS-true-tob",
    "NET-CS-true-tob",
]


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Boxplots from experiment 3 results CSV. "
        )
    )

    ap.add_argument(
        "--csv",
        required=True,
        help="Combined CSV, e.g. results-combined-exclude.csv",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Output directory for PDF figures",
    )
    ap.add_argument(
        "--taxa",
        nargs="*",
        type=int,
        default=None,
        help="Optional taxa subset, e.g. --taxa 50 100 150 200",
    )
    ap.add_argument(
        "--show-points",
        action="store_true",
        help="Overlay faint jittered replicate points",
    )

    return ap.parse_args()


def parse_taxa(value):
    if pd.isna(value):
        return np.nan
    m = re.search(r"(\d+)", str(value))
    return int(m.group(1)) if m else np.nan


def canonicalize_method(value):
    method = str(value).strip()

    mapping = {
        "ASTRAL": "ASTRAL",
        "CAMUS": "CAMUS",

        # True-ToB results
        "NANQUPLUS": "NANQUPLUS-true-tob",
        "NANQUPLUS-true-tob": "NANQUPLUS-true-tob",
        "NET-CS-true-tob": "NET-CS-true-tob",
        "NET-CS": "NET-CS-true-tob",

        # Estimated-ToB results
        "NANQUPLUS-TOB-QMC-3f1a": "NANQUPLUS-estimated-tob",
        "NET-CS-TOB-QMC-3f1a": "NET-CS-estimated-tob",

    }

    return mapping.get(method, method)


def get_method_family(method):
    if method == "ASTRAL":
        return "ASTRAL"
    if method == "CAMUS":
        return "CAMUS"
    if method.startswith("NANQUPLUS"):
        return "NANQUPLUS"
    if method.startswith("NET-CS"):
        return "NET-CS"
    return method


def get_source_variant(method):
    if method in {"ASTRAL", "CAMUS"}:
        return None
    if method.endswith("-estimated-tob"):
        return "estimated-tob"
    if method.endswith("-true-tob"):
        return "true-tob"
    return None


def get_method_color(method):
    if method in METHOD_COLORS:
        return METHOD_COLORS[method]

    family = get_method_family(method)
    if family in METHOD_FAMILY_COLORS:
        return METHOD_FAMILY_COLORS[family]

    return "black"


def prepare_dataframe(path):
    print("[INFO] Reading:", path)

    df = pd.read_csv(path)

    required = {
        "model",
        "replicate",
        "method",
        "normalized-hwcd",
        "leaves_hybrid_fn_rate",
        "leaves_hybrid_fp_rate",
    }

    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "Missing required column(s): " + ", ".join(sorted(missing))
        )

    df["taxa_num"] = df["model"].apply(parse_taxa)
    df["method"] = df["method"].apply(canonicalize_method)

    numeric_cols = [
        "normalized-hwcd",
        "leaves_hybrid_fn_rate",
        "leaves_hybrid_fp_rate",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["balanced_hybrid_error"] = (
        df["leaves_hybrid_fn_rate"]
        + df["leaves_hybrid_fp_rate"]
    ) / 2.0

    
    allowed = set(HWCD_METHOD_ORDER) | set(HYBRID_METHOD_ORDER)
    unexpected = sorted(set(df["method"]) - allowed)
    if unexpected:
        print(
            "[WARN] Ignoring unrecognized method(s): "
            + ", ".join(unexpected)
        )

    df = df[df["method"].isin(allowed)].copy()

    if df.empty:
        raise ValueError("No recognized plotting methods remain in the CSV.")

    return df


def sorted_unique(vals):
    return sorted(pd.unique(vals))


def method_block(method):
    if method == "CAMUS":
        return "camus"
    source = get_source_variant(method)
    return source if source is not None else method


def compute_method_offsets(method_order):
    if not method_order:
        return {}

    if len(method_order) == 1:
        return {method_order[0]: 0.0}

    raw = [0.0]
    previous_block = method_block(method_order[0])

    for method in method_order[1:]:
        current_block = method_block(method)
        step = 1.0 if current_block == previous_block else 1.55
        raw.append(raw[-1] + step)
        previous_block = current_block

    raw = np.asarray(raw, dtype=float)
    raw -= raw.mean()

    max_abs = np.max(np.abs(raw))
    if max_abs > 0:
        raw *= 0.47 / max_abs

    return dict(zip(method_order, raw))


def add_taxa_separators(ax, n_x):
    for i in range(n_x - 1):
        ax.axvline(
            x=(i + 0.5) * GROUP_GAP,
            color="0.95",
            linewidth=0.50,
            zorder=0,
        )


def style_axis(ax, ylabel):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_linewidth(0.7)
    ax.tick_params(axis="both", width=0.7, length=3)
    ax.set_ylabel(ylabel)


def has_missing_groups(
    plot_df,
    taxa_order,
    method_order,
    metric_col,
):
    for taxa in taxa_order:
        for method in method_order:
            vals = plot_df.loc[
                (plot_df["taxa_num"] == taxa)
                & (plot_df["method"] == method),
                metric_col,
            ].dropna()

            if len(vals) == 0:
                return True

    return False


def add_missing_group_marks(
    ax,
    plot_df,
    taxa_order,
    method_order,
    metric_col,
    y_fraction=0.035,
):
    offsets = compute_method_offsets(method_order)
    y0, y1 = ax.get_ylim()
    y = y0 + y_fraction * (y1 - y0)

    for i, taxa in enumerate(taxa_order):
        center = i * GROUP_GAP

        for method in method_order:
            vals = plot_df.loc[
                (plot_df["taxa_num"] == taxa)
                & (plot_df["method"] == method),
                metric_col,
            ].dropna()

            if len(vals) == 0:
                ax.text(
                    center + offsets[method],
                    y,
                    MISSING_GROUP_MARKER,
                    ha="center",
                    va="center",
                    fontsize=5.6,
                    color="0.45",
                    zorder=5,
                )


def grouped_boxplot(
    ax,
    data,
    metric,
    taxa_order,
    method_order,
    ylabel,
    show_points=False,
):
    sub = data[["taxa_num", "method", metric]].dropna()

    offsets = compute_method_offsets(method_order)
    box_w = min(0.68 / max(len(method_order), 1), 0.14)

    positions = []
    box_data = []
    methods_for_box = []

    for i, taxa in enumerate(taxa_order):
        center = i * GROUP_GAP

        for method in method_order:
            vals = sub.loc[
                (sub["taxa_num"] == taxa)
                & (sub["method"] == method),
                metric,
            ].dropna().to_numpy(dtype=float)

            if len(vals) == 0:
                continue

            positions.append(center + offsets[method])
            box_data.append(vals)
            methods_for_box.append(method)

    if box_data:
        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=box_w,
            patch_artist=True,
            showfliers=False,
            manage_ticks=False,
            medianprops={"color": "black", "linewidth": 1.0},
            boxprops={"linewidth": 0.75},
            whiskerprops={"linewidth": 0.75, "color": "0.35"},
            capprops={"linewidth": 0.75, "color": "0.35"},
        )

        for box, method in zip(bp["boxes"], methods_for_box):
            box.set_facecolor(get_method_color(method))
            box.set_edgecolor("0.25")
            box.set_linewidth(0.75)
            box.set_alpha(0.82)

        if show_points:
            rng = np.random.default_rng(12345)
            for pos, vals, method in zip(
                positions,
                box_data,
                methods_for_box,
            ):
                jitter = rng.normal(
                    0,
                    box_w * 0.10,
                    size=len(vals),
                )
                ax.scatter(
                    pos + jitter,
                    vals,
                    s=3,
                    alpha=0.12,
                    color=get_method_color(method),
                    edgecolors="none",
                    zorder=2,
                )

    ax.set_xticks(
        [i * GROUP_GAP for i in range(len(taxa_order))]
    )
    ax.set_xticklabels(
        [str(t) for t in taxa_order]
    )
    ax.set_xlabel("Number of taxa")

    style_axis(ax, ylabel)
    add_taxa_separators(ax, len(taxa_order))


def add_two_level_legend_to_axis(
    ax,
    present_methods,
    show_missing_note=False,
):
    from matplotlib.patches import Patch

    ax.set_axis_off()
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)

    # Preserve first appearance in hue order.
    present_families = []
    for method in present_methods:
        family = get_method_family(method)
        if family not in present_families:
            present_families.append(family)

    family_handles = [
        Patch(
            facecolor=METHOD_FAMILY_COLORS[family],
            edgecolor="0.25",
            linewidth=0.70,
            label=FAMILY_LEGEND_LABELS[family],
        )
        for family in present_families
    ]

    ax.text(
        0.5,
        0.90,
        "Reconstruction method",
        ha="center",
        va="center",
        fontsize=6.3,
        color="0.15",
    )

    if family_handles:
        family_legend = ax.legend(
            handles=family_handles,
            labels=[
                FAMILY_LEGEND_LABELS[f]
                for f in present_families
            ],
            loc="center",
            bbox_to_anchor=(0.5, 0.67),
            ncol=len(family_handles),
            frameon=False,
            borderaxespad=0.0,
            handlelength=1.25,
            handletextpad=0.45,
            columnspacing=1.35,
            fontsize=6.2,
        )
        ax.add_artist(family_legend)

    source_order = []
    for method in present_methods:
        source = get_source_variant(method)
        if source is not None and source not in source_order:
            source_order.append(source)

    source_handles = [
        Patch(
            facecolor=SOURCE_LEGEND_COLORS[source],
            edgecolor="0.25",
            linewidth=0.70,
            label=SOURCE_LEGEND_LABELS[source],
        )
        for source in source_order
    ]

    if source_handles:
        ax.text(
            0.5,
            0.40,
            "ToB source / variant",
            ha="center",
            va="center",
            fontsize=6.3,
            color="0.15",
        )

        ax.legend(
            handles=source_handles,
            labels=[
                SOURCE_LEGEND_LABELS[s]
                for s in source_order
            ],
            loc="center",
            bbox_to_anchor=(0.5, 0.16),
            ncol=len(source_handles),
            frameon=False,
            borderaxespad=0.0,
            handlelength=1.25,
            handletextpad=0.45,
            columnspacing=1.25,
            fontsize=6.2,
        )

    if show_missing_note:
        ax.text(
            0.985,
            0.16,
            f"{MISSING_GROUP_MARKER} = no valid result",
            ha="right",
            va="center",
            fontsize=5.9,
            color="0.50",
        )




def print_camus_vs_astral_hwcd(df, atol=1e-12):
    """
    Compare CAMUS vs ASTRAL on matched replicates using normalized HWCD.

    Lower HWCD is better.

    Prints better / tie / worse counts by taxa size and overall.
    Only replicate pairs where BOTH CAMUS and ASTRAL have a valid HWCD
    value are included.
    """
    metric = "normalized-hwcd"

    sub = df.loc[
        df["method"].isin(["CAMUS", "ASTRAL"]),
        ["taxa_num", "replicate", "method", metric],
    ].copy()

    sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
    sub = sub.dropna(subset=["taxa_num", "replicate", metric])

    if sub.empty:
        print()
        print("[CAMUS vs ASTRAL HWCD]")
        print("No matched CAMUS/ASTRAL HWCD data available.")
        return

    # If duplicate rows occur for the same taxa/replicate/method, use the first
    # after warning. In the expected input there should be one row per method.
    dup_mask = sub.duplicated(
        subset=["taxa_num", "replicate", "method"],
        keep=False,
    )
    if dup_mask.any():
        n_dup = int(dup_mask.sum())
        print(
            f"[WARN] CAMUS-vs-ASTRAL comparison found {n_dup} duplicate "
            "taxa/replicate/method rows; using the first occurrence."
        )
        sub = sub.drop_duplicates(
            subset=["taxa_num", "replicate", "method"],
            keep="first",
        )

    wide = sub.pivot(
        index=["taxa_num", "replicate"],
        columns="method",
        values=metric,
    )

    if "CAMUS" not in wide.columns or "ASTRAL" not in wide.columns:
        print()
        print("[CAMUS vs ASTRAL HWCD]")
        print("No matched CAMUS/ASTRAL HWCD data available.")
        return

    wide = wide.dropna(subset=["CAMUS", "ASTRAL"]).reset_index()

    if wide.empty:
        print()
        print("[CAMUS vs ASTRAL HWCD]")
        print("No replicates have valid HWCD for both CAMUS and ASTRAL.")
        return

    diff = wide["CAMUS"] - wide["ASTRAL"]

    wide["comparison"] = np.where(
        np.isclose(diff, 0.0, atol=atol, rtol=0.0),
        "tie",
        np.where(diff < 0.0, "better", "worse"),
    )

    print()
    print("[CAMUS vs ASTRAL HWCD]")
    print("Lower normalized HWCD is better.")
    print("Only matched replicates with valid results for both methods are counted.")

    summary_rows = []

    for taxa in sorted(pd.unique(wide["taxa_num"])):
        group = wide.loc[wide["taxa_num"] == taxa]

        better = int((group["comparison"] == "better").sum())
        tie = int((group["comparison"] == "tie").sum())
        worse = int((group["comparison"] == "worse").sum())

        summary_rows.append({
            "taxa": int(taxa),
            "matched_replicates": int(len(group)),
            "CAMUS_better": better,
            "tie": tie,
            "CAMUS_worse": worse,
        })

    overall_better = int((wide["comparison"] == "better").sum())
    overall_tie = int((wide["comparison"] == "tie").sum())
    overall_worse = int((wide["comparison"] == "worse").sum())

    summary_rows.append({
        "taxa": "ALL",
        "matched_replicates": int(len(wide)),
        "CAMUS_better": overall_better,
        "tie": overall_tie,
        "CAMUS_worse": overall_worse,
    })

    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))


def print_box_medians(
    df,
    metric,
    taxa_order,
    method_order,
    metric_label=None,
):
    """
    Print the median value used by every taxa/method box.

    Missing groups are printed as NR.
    """
    if metric_label is None:
        metric_label = metric

    print()
    print(f"[MEDIANS] {metric_label}")

    rows = []

    for taxa in taxa_order:
        for method in method_order:
            vals = df.loc[
                (df["taxa_num"] == taxa)
                & (df["method"] == method),
                metric,
            ].dropna()

            if len(vals) == 0:
                median = np.nan
                median_text = "NR"
            else:
                median = float(vals.median())
                median_text = f"{median:.6f}"

            rows.append({
                "taxa": int(taxa),
                "method": method,
                "n": int(len(vals)),
                "median": median_text,
            })

    median_df = pd.DataFrame(rows)

    print(
        median_df.to_string(
            index=False,
        )
    )


def save_metric_plot(
    df,
    metric,
    ylabel,
    out_pdf,
    y_limits,
    y_ticks,
    show_points=False,
):
    taxa_order = sorted_unique(
        df["taxa_num"].dropna()
    )

    existing = set(df["method"])

    if metric == "normalized-hwcd":
        base_method_order = HWCD_METHOD_ORDER
    else:
        base_method_order = HYBRID_METHOD_ORDER

    method_order = [
        m for m in base_method_order
        if m in existing
    ]

    # Print the exact median corresponding to every plotted box.
    print_box_medians(
        df=df,
        metric=metric,
        taxa_order=taxa_order,
        method_order=method_order,
        metric_label=ylabel,
    )

    fig = plt.figure(figsize=(7.5, 3.55))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.0, 0.11, 0.34],
        hspace=0.03,
    )

    ax = fig.add_subplot(gs[0])

    ax_spacer = fig.add_subplot(gs[1])
    ax_spacer.set_axis_off()

    ax_legend = fig.add_subplot(gs[2])

    grouped_boxplot(
        ax=ax,
        data=df,
        metric=metric,
        taxa_order=taxa_order,
        method_order=method_order,
        ylabel=ylabel,
        show_points=show_points,
    )

    # Hybrid metrics: retain the current slight negative margin so zero-valued
    # boxes/medians remain visually distinct from the x-axis.
    if metric in {
        "leaves_hybrid_fn_rate",
        "leaves_hybrid_fp_rate",
        "balanced_hybrid_error",
    }:
        lower = min(y_limits[0], -0.05)
        ax.set_ylim(lower, y_limits[1])
    else:
        ax.set_ylim(*y_limits)

    ax.set_yticks(y_ticks)

    if metric in {
        "leaves_hybrid_fn_rate",
        "leaves_hybrid_fp_rate",
        "balanced_hybrid_error",
    } and y_limits[1] >= 0.5:
        ax.axhline(
            y=0.5,
            color="0.55",
            linestyle="--",
            linewidth=0.7,
            alpha=0.85,
            zorder=1,
        )
        ax.text(
            0.995,
            0.505,
            "0.5",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=6,
            color="0.45",
        )

    missing = has_missing_groups(
        df,
        taxa_order,
        method_order,
        metric,
    )

    if missing:
        add_missing_group_marks(
            ax,
            df,
            taxa_order,
            method_order,
            metric,
            y_fraction=0.035,
        )

    add_two_level_legend_to_axis(
        ax_legend,
        method_order,
        show_missing_note=missing,
    )

    fig.subplots_adjust(
        left=0.085,
        right=0.99,
        top=0.975,
        bottom=0.035,
    )

    os.makedirs(
        os.path.dirname(out_pdf) or ".",
        exist_ok=True,
    )

    fig.savefig(
        out_pdf,
        bbox_inches="tight",
        dpi=600,
    )
    plt.close(fig)

    print("[INFO] Saved:", out_pdf)


def print_summary(df):
    print("\n[INFO] Rows by taxa and method:")
    counts = (
        df.groupby(["taxa_num", "method"], observed=True)
        .size()
        .rename("rows")
        .reset_index()
    )
    print(counts.to_string(index=False))


def main():
    args = parse_args()

    df = prepare_dataframe(args.csv)

    if args.taxa is not None:
        df = df[
            df["taxa_num"].isin(args.taxa)
        ].copy()

    if df.empty:
        raise ValueError("No rows remain after filtering.")

    print_summary(df)

    # Replicate-level matched comparison for normalized HWCD.
    print_camus_vs_astral_hwcd(df)

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    plots = [
        (
            "normalized-hwcd",
            "Normalized Hardwired Cluster Distance (HWCD)",
            "normalized_hwcd.pdf",
            (0.0, 0.25),
            np.arange(0.0, 0.251, 0.05),
        ),
        (
            "leaves_hybrid_fn_rate",
            "Hybrid FNR",
            "hybrid_fnr.pdf",
            (0.0, 1.0),
            np.arange(0.0, 1.01, 0.2),
        ),
        (
            "leaves_hybrid_fp_rate",
            "Hybrid FPR",
            "hybrid_fpr.pdf",
            (0.0, 1.0),
            np.arange(0.0, 1.01, 0.2),
        ),
        (
            "balanced_hybrid_error",
            "Balanced hybrid leaves error",
            "balanced_hybrid_error.pdf",
            (0.0, 1.0),
            np.arange(0.0, 1.01, 0.2),
        ),
    ]

    for metric, ylabel, filename, y_limits, y_ticks in plots:
        save_metric_plot(
            df=df,
            metric=metric,
            ylabel=ylabel,
            out_pdf=os.path.join(outdir, filename),
            y_limits=y_limits,
            y_ticks=y_ticks,
            show_points=args.show_points,
        )


if __name__ == "__main__":
    main()