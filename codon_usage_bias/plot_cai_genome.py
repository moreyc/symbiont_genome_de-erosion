import argparse
import csv
import math

import matplotlib.pyplot as plt
import pandas as pd

# Read the coreCAI tsv
def read_results(tsv):
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
                "interesting": row["interesting"] == "True"
            })
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Plot genome-wide CAI from the coreCAI TSV."
    )
    parser.add_argument(
        "-i", "--input",
        default="coreCAI.tsv",
        help="coreCAI TSV from core_cai.py"
    )
    parser.add_argument(
        "--plot",
        default="coreCAI_genome.png",
        help="Genome-wide CAI plot"
    )
    args = parser.parse_args()

    print("Reading results...")
    results = read_results(args.input)
    print(
        f"  Genes: {len(results)}"
    )

    # Genome-wide CAI plot
    print(
        "\nCreating genome-wide CAI plot..."
    )
    contigs = {}

    for r in results:
        if math.isfinite(r["CAI"]):
            contigs.setdefault(
                r["contig"],
                []
            ).append(r)

    fig, ax = plt.subplots(
        figsize=(14, 5)
    )

    if len(contigs) == 1:
        contig = list(contigs.keys())[0]
        data = sorted(
            contigs[contig],
            key=lambda r: r["position"]
        )
        normal = [
            r for r in data
            if not r["interesting"]
        ]
        interesting = [
            r for r in data
            if r["interesting"]
        ]

        rolling_window = 50
        y_rolling = pd.Series([r["CAI"] for r in data]).rolling(window=rolling_window, center=True).mean()

        if interesting:
            print("interesting")
            ax.scatter(
                [r["position"] / 1e6 for r in normal],
                [r["CAI"] for r in normal],
                s=60,
                alpha=0.6,
                label="All other genes",
                color="lightgray"
            )
            ax.scatter(
                [r["position"] / 1e6 for r in interesting],
                [r["CAI"] for r in interesting],
                s=60,
                alpha=0.6,
                label="ORFan genes",
                color="magenta"
            )
        else:
            ax.scatter(
                        [r["position"] / 1e6 for r in data],
                        [r["CAI"] for r in data],
                        alpha=0.6,
                        s=60,
                        color="lightgray"
                    )

        ax.plot(
                [r["position"] / 1e6 for r in data],
                y_rolling,
                color="black",
                linewidth=2,
                label=f"{rolling_window}-Point Rolling Avg",
            )

        # ax.set_xlabel(
        #     "Genome position (Mb)",
        #     fontsize=20
        # )

    else:
        offset = 0
        centers = []
        labels = []

        for contig, data in contigs.items():
            data = sorted(
                data,
                key=lambda r: r["position"]
            )
            length = max(
                r["end"] for r in data
            )
            normal = [
                r for r in data
                if not r["interesting"]
            ]
            interesting = [
                r for r in data
                if r["interesting"]
            ]
            ax.scatter(
                [
                    (r["position"] + offset) / 1e6
                    for r in normal
                ],
                [
                    r["CAI"]
                    for r in normal
                ],
                s=10,
                alpha=0.6
            )
            if interesting:
                ax.scatter(
                    [
                        (r["position"] + offset) / 1e6
                        for r in interesting
                    ],
                    [
                        r["CAI"]
                        for r in interesting
                    ],
                    s=10,
                    alpha=0.6,
                    label="Putative HGT region CDS"
                )
            centers.append(
                (offset + length / 2) / 1e6
            )
            labels.append(contig)
            offset += length

        ax.set_xticks(centers)
        ax.set_xticklabels(
            labels,
            rotation=45,
            ha="right",
            fontsize=20
        )
        ax.set_xlabel("Contig")

    # ax.set_ylabel(
    #     "CAI relative to conserved core",
    #     fontsize=20
    # )
    # ax.set_title(
    #     "Genome-wide codon adaptation index"
    # )
    ax.grid(
        axis="y",
        alpha=0.2
    )
    # if args.interesting_gff:
    #     ax.legend()

    plt.tight_layout()
    plt.savefig(
        args.plot,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    print(
        f"Saved: {args.plot}"
    )


if __name__ == "__main__":
    main()