#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ============================================================
# Default experiment layout
# ============================================================

DEFAULT_TAXA = [50, 100, 200]
DEFAULT_GT_COUNTS = [100, 500, 1000]
DEFAULT_ILS = [0.25, 0.5, 1.0, 2.0]


def extract_num_gene_trees(method):
    """
    Convert:

        NET-CS-100-GT
        NET-CS-500-GT
        NET-CS-1000-GT

    into:

        100
        500
        1000

    Matching is case-insensitive.
    """
    if pd.isna(method):
        return np.nan

    s = str(method).strip()

    m = re.fullmatch(
        r"NET-CS-(100|500|1000)-GT",
        s,
        flags=re.IGNORECASE,
    )

    if m is None:
        return np.nan

    return int(m.group(1))


def ils_label(value):
    value = float(value)

    if value.is_integer():
        return f"{value:.1f}x"

    return f"{value:g}x"


def get_metric_info(metric):
    """
    Return:
        source column
        y-axis label
        whether conversion is required
    """

    metric = metric.lower()

    if metric == "kendall":
        return (
            "avg_kendall_norm",
            "Average Kendall-tau distance",
            None,
        )

    if metric == "hybrid":
        return (
            "hybrid_error_rate",
            "Average hybrid error rate",
            None,
        )

    if metric == "runtime":
        return (
            "runtime(seconds)",
            "Runtime (seconds)",
            None,
        )

    if metric == "rss":
        return (
            "peak_RSS",
            "Peak RSS (MB)",
            "kb_to_mb",
        )

    raise ValueError(
        f"Unknown metric: {metric}. "
        "Choose from: kendall, hybrid, runtime, rss"
    )


def plot_gt_boxplots(
    df,
    output,
    metric="kendall",
    taxa_values=None,
    gt_counts=None,
    ils_values=None,
    show_fliers=False,
):
    if taxa_values is None:
        taxa_values = DEFAULT_TAXA

    if gt_counts is None:
        gt_counts = DEFAULT_GT_COUNTS

    if ils_values is None:
        ils_values = DEFAULT_ILS

    metric_col, ylabel, conversion = get_metric_info(metric)

    # ========================================================
    # Prepare numeric GT column
    # ========================================================

    df = df.copy()

    df["num_gene_trees"] = (
        df["method"]
        .apply(extract_num_gene_trees)
    )

    df = df[
        df["num_gene_trees"].notna()
    ].copy()

    df["num_gene_trees"] = (
        df["num_gene_trees"].astype(int)
    )

    # ========================================================
    # Numeric conversions
    # ========================================================

    for col in [
        "taxa(round)",
        "ils_level",
        metric_col,
    ]:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce",
        )

    # RSS is stored as kB.
    if conversion == "kb_to_mb":
        df["_plot_metric"] = (
            df[metric_col] / 1024.0
        )
    else:
        df["_plot_metric"] = df[metric_col]

    # ========================================================
    # Keep requested experiment conditions
    # ========================================================

    df = df[
        df["taxa(round)"].isin(taxa_values)
        & df["num_gene_trees"].isin(gt_counts)
        & df["ils_level"].isin(ils_values)
    ].copy()

    df = df.dropna(
        subset=[
            "taxa(round)",
            "ils_level",
            "num_gene_trees",
            "_plot_metric",
        ]
    )

    if df.empty:
        raise ValueError(
            "No usable rows remain after filtering."
        )

    # ========================================================
    # Diagnostics
    # ========================================================

    print("Rows used:")
    print(
        df.groupby(
            [
                "taxa(round)",
                "num_gene_trees",
                "ils_level",
            ]
        )
        .size()
        .to_string()
    )

    # ========================================================
    # Figure
    # ========================================================

    n_panels = len(taxa_values)

    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(15.2, 3.7),
        squeeze=False,
    )

    axes = axes[0]

    # Matplotlib default qualitative palette, close to the
    # screenshot: blue, orange, green, red.
    colors = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
    ]

    # Center of each GT group.
    x_centers = np.arange(
        len(gt_counts),
        dtype=float,
    )

    # Width of each individual box.
    box_width = 0.16

    # Horizontal offsets for ILS levels.
    #
    # For four levels:
    #   -0.27, -0.09, +0.09, +0.27
    #
    spacing = 0.18

    offsets = (
        np.arange(len(ils_values))
        - (len(ils_values) - 1) / 2.0
    ) * spacing

    panel_letters = [
        chr(ord("A") + i)
        for i in range(n_panels)
    ]

    # ========================================================
    # Draw each taxa panel
    # ========================================================

    for panel_idx, (ax, taxa) in enumerate(
        zip(axes, taxa_values)
    ):

        taxa_df = df[
            df["taxa(round)"] == taxa
        ]

        # ----------------------------------------------------
        # One boxplot call per ILS level
        # ----------------------------------------------------

        for ils_idx, ils in enumerate(ils_values):

            box_data = []
            positions = []

            for gt_idx, gt in enumerate(gt_counts):

                values = (
                    taxa_df[
                        (taxa_df["num_gene_trees"] == gt)
                        & (taxa_df["ils_level"] == ils)
                    ]["_plot_metric"]
                    .dropna()
                    .to_numpy()
                )

                # matplotlib.boxplot cannot use an empty array.
                # We simply skip missing combinations.
                if len(values) == 0:
                    continue

                box_data.append(values)

                positions.append(
                    x_centers[gt_idx]
                    + offsets[ils_idx]
                )

            if box_data:

                bp = ax.boxplot(
                    box_data,
                    positions=positions,
                    widths=box_width,
                    patch_artist=True,
                    showfliers=show_fliers,
                    manage_ticks=False,
                    medianprops={
                        "linewidth": 1.25,
                    },
                    whiskerprops={
                        "linewidth": 1.0,
                    },
                    capprops={
                        "linewidth": 1.0,
                    },
                    boxprops={
                        "linewidth": 1.0,
                    },
                )

                for patch in bp["boxes"]:
                    patch.set_facecolor(
                        colors[ils_idx]
                    )
                    patch.set_alpha(0.55)

        # ----------------------------------------------------
        # Axes
        # ----------------------------------------------------

        ax.set_xticks(x_centers)

        ax.set_xticklabels(
            [str(x) for x in gt_counts]
        )

        ax.set_title(
            f"({panel_letters[panel_idx]}) "
            f"{taxa} taxon",
            fontsize=12,
            pad=10,
        )

        ax.grid(
            axis="y",
            alpha=0.20,
            linewidth=0.8,
        )

        ax.set_axisbelow(True)

        # Small horizontal padding.
        ax.set_xlim(
            -0.5,
            len(gt_counts) - 0.5,
        )

        ax.tick_params(
            axis="both",
            labelsize=10,
        )

        # Only first panel needs ylabel.
        if panel_idx == 0:
            ax.set_ylabel(
                ylabel,
                fontsize=11,
            )

    # ========================================================
    # Shared legend
    # ========================================================

    legend_handles = [
        Patch(
            facecolor=colors[i],
            alpha=0.55,
            label=ils_label(ils),
        )
        for i, ils in enumerate(ils_values)
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=len(ils_values),
        frameon=False,
        fontsize=10,
        handlelength=2.8,
        columnspacing=2.5,
    )

    # ========================================================
    # Layout
    # ========================================================

    fig.subplots_adjust(
        left=0.075,
        right=0.985,
        top=0.88,
        bottom=0.22,
        wspace=0.18,
    )

    # ========================================================
    # Write PDF
    # ========================================================

    os.makedirs(
        os.path.dirname(
            os.path.abspath(output)
        ),
        exist_ok=True,
    )

    fig.savefig(
        output,
        bbox_inches="tight",
    )

    plt.close(fig)

    print()
    print(f"Plot written to:")
    print(f"  {output}")


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Create grouped boxplots for the NET-CS "
            "100/500/1000 gene-tree experiment."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input GT CSV",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output PDF",
    )

    parser.add_argument(
        "--metric",
        choices=[
            "kendall",
            "hybrid",
            "runtime",
            "rss",
        ],
        default="kendall",
        help=(
            "Metric to plot. Default: kendall"
        ),
    )

    parser.add_argument(
        "--show-fliers",
        action="store_true",
        help="Show boxplot outliers",
    )

    args = parser.parse_args()

    # ========================================================
    # Read CSV
    # ========================================================

    input_path = os.path.abspath(
        args.input
    )

    df = pd.read_csv(input_path)

    # In case column names contain escaped underscores.
    df.columns = [
        str(col).replace("\\_", "_")
        for col in df.columns
    ]

    # ========================================================
    # Validate common columns
    # ========================================================

    required = [
        "taxa(round)",
        "ils_level",
        "method",
    ]

    missing = [
        col
        for col in required
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
        )

    # ========================================================
    # Default filename
    # ========================================================

    if args.output is None:

        input_dir = os.path.dirname(
            input_path
        )

        input_stem = os.path.splitext(
            os.path.basename(input_path)
        )[0]

        args.output = os.path.join(
            input_dir,
            f"{input_stem}_{args.metric}_boxplot.pdf",
        )

    plot_gt_boxplots(
        df=df,
        output=args.output,
        metric=args.metric,
        taxa_values=[
            50,
            100,
            200,
        ],
        gt_counts=[
            100,
            500,
            1000,
        ],
        ils_values=[
            0.25,
            0.5,
            1.0,
            2.0,
        ],
        show_fliers=args.show_fliers,
    )


if __name__ == "__main__":
    main()