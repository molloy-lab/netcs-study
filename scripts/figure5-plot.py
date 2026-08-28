#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def alpha_sort_key(x):
    try:
        return float(x)
    except Exception:
        return str(x)


def grouped_boxplot(ax, df, value_col, ylabel, method_order, alpha_order, method_colors):
    methods_present = [m for m in method_order if m in set(df["method"])]

    if not methods_present or not alpha_order:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_axis_off()
        return

    group_centers = np.arange(len(alpha_order)) * 1.5
    n_methods = len(methods_present)
    box_width = min(0.7 / max(n_methods, 1), 0.3)

    if n_methods == 1:
        offsets = [0.0]
    else:
        offsets = np.linspace(-0.22, 0.22, n_methods)

    legend_handles = []

    for j, method in enumerate(methods_present):
        data_list = []
        positions = []

        for i, alpha in enumerate(alpha_order):
            vals = df.loc[
                (df["method"] == method) & (df["alpha_str"] == alpha),
                value_col
            ].dropna().values

            if len(vals) > 0:
                data_list.append(vals)
                positions.append(group_centers[i] + offsets[j])

        if not data_list:
            continue

        color = method_colors.get(method, "gray")

        bp = ax.boxplot(
            data_list,
            positions=positions,
            widths=box_width,
            patch_artist=True,
            manage_ticks=False
        )

        for box in bp["boxes"]:
            box.set_facecolor(color)
            box.set_alpha(0.8)
            box.set_linewidth(1.1)

        for whisker in bp["whiskers"]:
            whisker.set_linewidth(1.1)

        for cap in bp["caps"]:
            cap.set_linewidth(1.1)

        for median in bp["medians"]:
            median.set_linewidth(1.5)

        legend_handles.append(Patch(facecolor=color, edgecolor="black", label=method, alpha=0.8))

    ax.set_xticks(group_centers)
    ax.set_xticklabels(alpha_order, rotation=45, ha="right")
    ax.set_xlabel("alpha")
    ax.set_ylabel(ylabel)
    ax.legend(handles=legend_handles)
    ax.grid(axis="y", linestyle="--", alpha=0.5)


def main():
    parser = argparse.ArgumentParser(
        description="figure 5 plot."
    )
    parser.add_argument("--csv", required=True, help="Input CSV file")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)

    required_cols = ["alpha", "method", "(fnr+fpr)/2"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for col in ["(fnr+fpr)/2"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["alpha_str"] = df["alpha"].astype(str)

    method_order = ["TOB-QMC-3f1a", "TOB-QMC-default"]
    method_colors = {
        "TOB-QMC-3f1a": "#4C78A8",
        "TOB-QMC-default": "#F58518",
    }

    alpha_order = sorted(df["alpha_str"].dropna().unique(), key=alpha_sort_key)


    
    fig, ax = plt.subplots(figsize=(10, 5))
    grouped_boxplot(
        ax=ax,
        df=df,
        value_col="(fnr+fpr)/2",
        ylabel="Arithmetic mean of FNR and FPR",
        method_order=method_order,
        alpha_order=alpha_order,
        method_colors=method_colors,
    )
    fig.tight_layout()
    rf_out = f"figure5-draft.pdf"
    fig.savefig(rf_out)
    plt.close(fig)

    print(f"Saved: {rf_out}")


if __name__ == "__main__":
    main()