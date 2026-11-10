import argparse
import csv
import math

import matplotlib.pyplot as plt
import numpy as np


def parse_gff(gff_file):
    cds = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            (
                contig,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes
            ) = fields
            if feature != "CDS":
                continue
            attrs = {}
            for item in attributes.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            gene_id = (
                attrs.get("ID")
                or attrs.get("locus_tag")
                or attrs.get("gene")
                or attrs.get("Parent")
            )
            if gene_id is None:
                continue
            cds[gene_id] = {
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "product": attrs.get("product", "")
            }
    return cds

# Read the coreCAI tsv
def read_results(tsv, interesting_ids):
    def to_float(x):
        if x == "NA":
            return float("nan")
        return float(x)

    results = []
    with open(tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            results.append({
                "gene_id": row["gene_id"],
                "contig": row["contig"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "position": float(row["position"]),
                "CAI": to_float(row["coreCAI"]),
                "GC": to_float(row["GC"]),
                "GC3s": to_float(row["GC3s"]),
                "interesting":
                    row["gene_id"] in interesting_ids
            })
    return results


def add_regression(ax, x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if len(x) < 2:
        return
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(
        x.min(),
        x.max(),
        100
    )
    y_line = slope * x_line + intercept
    ax.plot(
        x_line,
        y_line,
        linestyle="--",
        linewidth=1.5,
        label="Linear regression",
        color="black"
    )


def make_scatter_plot(
    results,
    x_key,
    x_label,
    output,
    title,
    interesting_label="Putative HGT region CDS"
):
    normal = [
        r for r in results
        if not r["interesting"]
        and math.isfinite(r[x_key])
        and math.isfinite(r["CAI"])
    ]
    interesting = [
        r for r in results
        if r["interesting"]
        and math.isfinite(r[x_key])
        and math.isfinite(r["CAI"])
    ]

    fig, ax = plt.subplots(
        figsize=(6, 6)
    )

    # Normal genes
    ax.scatter(
        [r[x_key] for r in normal],
        [r["CAI"] for r in normal],
        s=12,
        alpha=0.4,
        label="All other genes",
        color="lightgray"
    )

    # Interesting genes
    if interesting:
        ax.scatter(
            [r[x_key] for r in interesting],
            [r["CAI"] for r in interesting],
            s=35,
            alpha=0.9,
            label=interesting_label,
            color="orange"
        )

    # Regression across ALL genes
    x = [
        r[x_key]
        for r in results
        if math.isfinite(r[x_key])
        and math.isfinite(r["CAI"])
    ]
    y = [
        r["CAI"]
        for r in results
        if math.isfinite(r[x_key])
        and math.isfinite(r["CAI"])
    ]
    add_regression(
        ax,
        x,
        y,
    )

    ax.set_xlabel(x_label, fontsize=18)
    ax.set_ylabel(
        "CAI relative to conserved core", fontsize=18
    )
    # ax.set_title(title)
    ax.grid(
        axis="both",
        alpha=0.2
    )
    # ax.legend(fontsize=18)

    plt.tight_layout()
    plt.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    print(f"Saved: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot CAI vs GC3s from the coreCAI TSV."
    )
    parser.add_argument(
        "-i", "--input",
        default="coreCAI.tsv",
        help="coreCAI TSV from core_cai.py"
    )
    parser.add_argument(
        "--gc3-plot",
        default="coreCAI_vs_GC3s.png",
        help="CAI vs GC3s plot"
    )
    parser.add_argument(
        "--gff",
        default=None,
        help="Genome GFF3 (required with --interesting-core)"
    )
    parser.add_argument(
        "--interesting-gff",
        default=None,
        help="Optional GFF3 containing interesting genes"
    )
    parser.add_argument(
        "--interesting-core",
        action="store_true",
        help="Flag the reference single-copy orthologs as interesting genes"
    )
    parser.add_argument(
        "--interesting-label",
        default=None,
        help="Legend label for interesting genes"
    )
    args = parser.parse_args()

    if args.interesting_core and not args.gff:
        parser.error("--interesting-core requires --gff")

    # Read interesting genes
    interesting_ids = set()

    if args.interesting_gff:
        print("\nReading interesting-gene GFF...")
        interesting_gff = parse_gff(args.interesting_gff)
        interesting_ids |= set(interesting_gff.keys())
        print(f"  Found {len(interesting_ids)} interesting genes.")

    if args.interesting_core:
        # core_ids = set(key.split("|")[1] for key in reference.keys())
        # non_core = set(genome.keys()) - core_ids
        # matched = core_ids & set(genome.keys())
        gff = parse_gff(args.gff)
        ribosomal = {
            gene_id for gene_id, info in gff.items()
            if "ribosomal" in info["product"].lower()
            }
        print(f"\nFound {len(ribosomal)} CDSs with 'ribosomal' in the product name.")
        interesting_ids |= ribosomal
        # print(
        #     f"\nFlagging single-copy orthologs as interesting: "
        #     f"{len(matched)} of {len(core_ids)} reference genes "
        #     "found in genome CDSs."
        # )

    if args.interesting_label:
        interesting_label = args.interesting_label
    elif args.interesting_core and not args.interesting_gff:
        interesting_label = "Single-copy orthologs"
    else:
        interesting_label = "Putative HGT region CDS"

    print("\nReading results...")
    results = read_results(args.input, interesting_ids)
    print(
        f"  Genes: {len(results)}"
    )

    # CAI vs GC3s
    print(
        "\nCreating CAI vs GC3s plot..."
    )
    make_scatter_plot(
        results,
        x_key="GC3s",
        x_label="GC3s",
        output=args.gc3_plot,
        title="CAI vs GC3s",
        interesting_label=interesting_label
    )


if __name__ == "__main__":
    main()