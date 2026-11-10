import argparse
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Test whether genes in a region sit above or below the "
            "genome-wide coreCAI vs GC3s regression line, using "
            "same-size sliding windows as the null."
        )
    )
    parser.add_argument("tsv", help="Output TSV from the CAI script")
    parser.add_argument(
        "region",
        help="contig:start-end in GFF coordinates (1-based, inclusive)"
    )
    args = parser.parse_args()

    contig, coords = args.region.rsplit(":", 1)
    rstart, rend = map(int, coords.replace(",", "").split("-"))
    size = rend - rstart + 1

    # Fit coreCAI ~ GC3s across all genes and compute residuals
    df = pd.read_csv(args.tsv, sep="\t", na_values="NA")
    df = df.dropna(subset=["coreCAI", "GC3s"]).copy()
    slope, intercept = np.polyfit(df["GC3s"], df["coreCAI"], 1)
    df["resid"] = df["coreCAI"] - (slope * df["GC3s"] + intercept)

    # Observed mean residual of genes whose midpoint is in the region
    in_region = (df["contig"] == contig) & df["position"].between(rstart, rend)
    if not in_region.any():
        sys.exit(f"No genes found in {args.region}")
    obs = df.loc[in_region, "resid"].mean()

    # Null: mean residual in same-size windows slid along every contig,
    # skipping windows that overlap the region itself
    step = int((df["end"] - df["start"] + 1).median())
    null = []

    for name, g in df.groupby("contig"):
        g = g.sort_values("position")
        pos = g["position"].to_numpy()
        csum = np.concatenate([[0], np.cumsum(g["resid"].to_numpy())])

        starts = np.arange(1, int(g["end"].max()) - size + 2, step)
        ends = starts + size - 1

        if name == contig:
            keep = (ends < rstart) | (starts > rend)
            starts, ends = starts[keep], ends[keep]

        lo = np.searchsorted(pos, starts, side="left")
        hi = np.searchsorted(pos, ends, side="right")
        n = hi - lo
        has_genes = n > 0
        null.append((csum[hi] - csum[lo])[has_genes] / n[has_genes])

    null = np.concatenate(null)

    if len(null) == 0:
        sys.exit("No null windows available for this region size")

    # Empirical p-values with a +1 correction
    p_above = (1 + np.sum(null >= obs)) / (len(null) + 1)
    p_below = (1 + np.sum(null <= obs)) / (len(null) + 1)
    p_two = (
        1 + np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean()))
    ) / (len(null) + 1)

    print(f"Region:               {args.region} ({size} bp)")
    print(f"Genes in region:      {int(in_region.sum())}")
    print(f"Mean gene length:     {int((df["end"] - df["start"] + 1).median())} bp")
    print(f"Regression:           CAI = {slope:.4f} * GC3s + {intercept:.4f}")
    print(f"Observed mean resid:  {obs:.5f}")
    print(f"Null windows:         {len(null)}")
    print(f"Null mean +/- SD:     {null.mean():.5f} +/- {null.std():.5f}")
    print(f"Z-score:              {(obs - null.mean()) / null.std():.2f}")
    print(f"P (above line):       {p_above:.4g}")
    print(f"P (below line):       {p_below:.4g}")
    print(f"P (two-sided):        {p_two:.4g}")


if __name__ == "__main__":
    main()