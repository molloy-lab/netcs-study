#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ============================================================
# Configuration
# ============================================================

TAXA_ORDER = [50, 100, 200]
ILS_ORDER = [0.25, 0.5, 1.0, 2.0]

METHOD_LABELS = {
    "NET-CS-TRUE-TOB": "Net-CS",
    "NANQUPLUS-TRUE-TOB": "NANUQ+",
}

METHOD_COLORS = {
    "NET-CS-TRUE-TOB": "#6ba6cf",
    "NANQUPLUS-TRUE-TOB": "#f5a45d",
}

PANEL_LABELS = {
    50: "(A) 50 taxon",
    100: "(B) 100 taxon",
    200: "(C) 200 taxon",
}

# Expected methods by taxa for complete-case plotting.
EXPECTED_METHODS_BY_TAXA = {
    50: ["NET-CS-TRUE-TOB", "NANQUPLUS-TRUE-TOB"],
    100: ["NET-CS-TRUE-TOB", "NANQUPLUS-TRUE-TOB"],
    200: ["NET-CS-TRUE-TOB"],
}


# ============================================================
# Helpers
# ============================================================

def normalize_column_name(name):
    name = str(name).strip()
    name = name.replace("\\_", "_")
    name = name.replace("**", "")
    return name.strip()


def read_results(path):
    df = pd.read_csv(path)

    df.columns = [normalize_column_name(c) for c in df.columns]

    required = [
        "taxa(round)",
        "network_id",
        "ils_level",
        "method",
        "avg_kendall_norm",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns: " + ", ".join(missing)
        )

    df["taxa(round)"] = pd.to_numeric(df["taxa(round)"], errors="coerce")
    df["ils_level"] = pd.to_numeric(df["ils_level"], errors="coerce")
    df["avg_kendall_norm"] = pd.to_numeric(df["avg_kendall_norm"], errors="coerce")

    df["network_id"] = df["network_id"].astype(str).str.strip()
    df["method"] = df["method"].astype(str).str.strip()

    df = df[
        df["taxa(round)"].isin(TAXA_ORDER)
        & df["ils_level"].isin(ILS_ORDER)
        & df["method"].isin(
            sorted({m for ms in EXPECTED_METHODS_BY_TAXA.values() for m in ms})
        )
        & np.isfinite(df["avg_kendall_norm"])
    ].copy()

    return df


def complete_case_filter(df):
    """
    For each (taxa, ils_level):
      - determine the expected methods for that taxa,
      - if any expected method is completely absent, drop the whole condition,
      - otherwise keep only network_id values shared by all expected methods.
    """
    pieces = []
    dropped_conditions = []

    for taxa in TAXA_ORDER:
        expected_methods = EXPECTED_METHODS_BY_TAXA[taxa]

        for ils in ILS_ORDER:
            cond = df[
                (df["taxa(round)"] == taxa)
                & np.isclose(df["ils_level"], ils)
            ].copy()

            if cond.empty:
                dropped_conditions.append((taxa, ils, "no rows"))
                continue

            methods_present = set(cond["method"].unique())

            if not set(expected_methods).issubset(methods_present):
                missing_methods = sorted(set(expected_methods) - methods_present)
                dropped_conditions.append(
                    (taxa, ils, f"missing required method(s): {', '.join(missing_methods)}")
                )
                continue

            id_sets = []
            for method in expected_methods:
                ids = set(
                    cond.loc[cond["method"] == method, "network_id"]
                )
                id_sets.append(ids)

            common_ids = set.intersection(*id_sets) if id_sets else set()

            if len(common_ids) == 0:
                dropped_conditions.append((taxa, ils, "no shared network_id across methods"))
                continue

            cond = cond[cond["network_id"].isin(common_ids)].copy()
            pieces.append(cond)

    if pieces:
        out = pd.concat(pieces, ignore_index=True)
    else:
        out = df.iloc[0:0].copy()

    return out, dropped_conditions


def print_diagnostics(df, dropped_conditions):
    print()
    print("Conditions used in plot")
    print("=" * 72)

    for taxa in TAXA_ORDER:
        print()
        print(f"{taxa} taxa")
        expected_methods = EXPECTED_METHODS_BY_TAXA[taxa]
        print(f"  expected methods: {', '.join(expected_methods)}")

        for ils in ILS_ORDER:
            cond = df[
                (df["taxa(round)"] == taxa)
                & np.isclose(df["ils_level"], ils)
            ]

            if cond.empty:
                print(f"  ILS {ils:g}x: DROPPED")
                continue

            shared_n = cond["network_id"].nunique()
            print(f"  ILS {ils:g}x: kept ({shared_n} shared replicate(s))")

            for method in expected_methods:
                n = cond.loc[cond["method"] == method, "network_id"].nunique()
                print(f"    {METHOD_LABELS[method]}: {n}")

    if dropped_conditions:
        print()
        print("Dropped conditions")
        print("=" * 72)
        for taxa, ils, reason in dropped_conditions:
            print(f"taxa={taxa}, ils={ils:g}: {reason}")


def get_offsets_for_methods(methods):
    """
    Center boxes depending on how many methods are plotted in that panel.
    """
    if len(methods) == 1:
        return {methods[0]: 0.0}
    if len(methods) == 2:
        return {
            methods[0]: -0.15,
            methods[1]: +0.15,
        }

    # fallback for >2 methods if ever needed
    span = np.linspace(-0.25, 0.25, len(methods))
    return {m: x for m, x in zip(methods, span)}


def make_plot(df, output, ymin=-0.005, ymax=0.112, show_fliers=False):
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.6))
    x = np.arange(len(ILS_ORDER))
    box_width = 0.18

    for ax, taxa in zip(axes, TAXA_ORDER):
        taxa_df = df[df["taxa(round)"] == taxa].copy()
        expected_methods = EXPECTED_METHODS_BY_TAXA[taxa]

        # only plot conditions that survived complete-case filtering
        valid_ils = []
        for ils in ILS_ORDER:
            cond = taxa_df[np.isclose(taxa_df["ils_level"], ils)]
            if not cond.empty:
                valid_ils.append(ils)

        offsets = get_offsets_for_methods(expected_methods)

        for method in expected_methods:
            method_df = taxa_df[taxa_df["method"] == method].copy()
            if method_df.empty:
                continue

            data = []
            positions = []

            for i, ils in enumerate(ILS_ORDER):
                cond_vals = method_df.loc[
                    np.isclose(method_df["ils_level"], ils),
                    "avg_kendall_norm"
                ].dropna().to_numpy()

                # if the whole condition was dropped, skip it
                if len(cond_vals) == 0:
                    continue

                data.append(cond_vals)
                positions.append(x[i] + offsets[method])

            if not data:
                continue

            bp = ax.boxplot(
                data,
                positions=positions,
                widths=box_width,
                patch_artist=True,
                showfliers=show_fliers,
                whis=1.5,
                medianprops={
                    "color": "#ff8c22",
                    "linewidth": 1.5,
                },
                boxprops={
                    "edgecolor": "#666666",
                    "linewidth": 0.8,
                },
                whiskerprops={
                    "color": "#444444",
                    "linewidth": 1.0,
                },
                capprops={
                    "color": "#444444",
                    "linewidth": 1.0,
                },
            )

            for box in bp["boxes"]:
                box.set_facecolor(METHOD_COLORS[method])
                box.set_alpha(0.9)

        ax.set_title(PANEL_LABELS[taxa], fontsize=22, pad=14)

        ax.set_xticks(x)
        ax.set_xticklabels(["0.25x", "0.5x", "1.0x", "2.0x"], fontsize=15)

        ax.set_xlim(-0.5, len(ILS_ORDER) - 0.5)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(np.arange(0.0, ymax + 1e-9, 0.02))
        ax.tick_params(axis="y", labelsize=15)

        ax.grid(axis="y", linestyle="-", linewidth=0.5, alpha=0.18)
        ax.grid(axis="x", visible=False)

        for spine in ax.spines.values():
            spine.set_linewidth(0.9)

    axes[0].set_ylabel("Average Kendall-tau Distance", fontsize=18)

    legend_handles = [
        Patch(
            facecolor=METHOD_COLORS["NET-CS-TRUE-TOB"],
            edgecolor="none",
            label="Net-CS",
        ),
        Patch(
            facecolor=METHOD_COLORS["NANQUPLUS-TRUE-TOB"],
            edgecolor="none",
            label="NANUQ+",
        ),
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=2,
        frameon=False,
        fontsize=18,
        handlelength=2.6,
        handleheight=1.0,
        columnspacing=5.0,
    )

    fig.subplots_adjust(
        left=0.065,
        right=0.992,
        top=0.86,
        bottom=0.26,
        wspace=0.18,
    )

    output = os.path.abspath(output)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    print()
    print(f"Wrote figure: {output}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot boxplots of avg_kendall_norm by taxa, ILS, and method "
            "using complete-case filtering by condition."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV file"
    )

    parser.add_argument(
        "-o", "--output",
        default="experiment2_kendall_boxplot_complete_case.pdf",
        help="Output figure file (default: experiment2_kendall_boxplot_complete_case.pdf)"
    )

    parser.add_argument(
        "--show-fliers",
        action="store_true",
        help="Show boxplot outliers"
    )

    parser.add_argument(
        "--ymin",
        type=float,
        default=-0.005,
        help="Y-axis minimum (default: -0.005)"
    )

    parser.add_argument(
        "--ymax",
        type=float,
        default=0.112,
        help="Y-axis maximum (default: 0.112)"
    )

    args = parser.parse_args()

    df = read_results(args.input)
    print(f"Rows after basic filtering: {len(df)}")

    df_cc, dropped_conditions = complete_case_filter(df)
    print(f"Rows after complete-case condition filtering: {len(df_cc)}")

    print_diagnostics(df_cc, dropped_conditions)

    make_plot(
        df_cc,
        output=args.output,
        ymin=args.ymin,
        ymax=args.ymax,
        show_fliers=args.show_fliers,
    )


if __name__ == "__main__":
    main()