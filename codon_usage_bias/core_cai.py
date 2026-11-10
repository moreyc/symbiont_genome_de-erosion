import argparse
import csv
import math

from Bio import SeqIO
from CAI import CAI, relative_adaptiveness


def read_fasta(fasta):
    return {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(fasta, "fasta")
    }

# Parse CDS features from gff
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


def calculate_gc(seq):
    seq = seq.upper()
    valid = [
        base for base in seq
        if base in "ACGT"
    ]
    if not valid:
        return float("nan")
    return (
        sum(base in "GC" for base in valid)
        / len(valid)
    )


def calculate_gc3s(seq):
    seq = seq.upper()

    # Genetic code
    codon_table = {
        "TTT": "F", "TTC": "F",
        "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I",
        "ATG": "M",
        "GTT": "V", "GTC": "V",
        "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S",
        "TCA": "S", "TCG": "S",
        "AGT": "S", "AGC": "S",
        "CCT": "P", "CCC": "P",
        "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T",
        "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y",
        "CAT": "H", "CAC": "H",
        "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N",
        "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C",
        "TGG": "W",
        "CGT": "R", "CGC": "R",
        "CGA": "R", "CGG": "R",
        "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G",
        "GGA": "G", "GGG": "G"
    }

    # Group synonymous codons
    aa_codons = {}
    for codon, aa in codon_table.items():
        aa_codons.setdefault(aa, []).append(codon)

    gc3 = 0
    total = 0
    usable_length = len(seq) - (len(seq) % 3)

    for i in range(0, usable_length, 3):
        codon = seq[i:i + 3]
        if codon not in codon_table:
            continue
        aa = codon_table[codon]
        # Only synonymous codon families
        if len(aa_codons[aa]) <= 1:
            continue
        third_base = codon[2]
        if third_base in "GC":
            gc3 += 1
        total += 1

    if total == 0:
        return float("nan")
    return gc3 / total


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate custom CAI, GC, and GC3s for all "
            "CDSs and write them to a TSV."
        )
    )
    parser.add_argument(
        "-g", "--genome",
        required=True,
        help="Full genome CDS FASTA (.ffn)"
    )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Focal-genome single-copy ortholog FASTA"
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Genome GFF3"
    )
    parser.add_argument(
        "--interesting-gff",
        default=None,
        help="Optional GFF3 containing interesting genes"
    )
    parser.add_argument(
        "-o", "--output",
        default="coreCAI.tsv",
        help="Output TSV"
    )
    parser.add_argument(
        "--interesting-core",
        action="store_true",
        help="Flag the reference single-copy orthologs as interesting genes"
    )
    args = parser.parse_args()


    # Read sequences
    print("Reading reference sequences...")
    reference = read_fasta(args.reference)
    print(
        f"  Reference genes: {len(reference)}"
    )
    print("Reading genome CDSs...")
    genome = read_fasta(args.genome)
    print(
        f"  Genome CDSs: {len(genome)}"
    )

    # Calculate reference codon weights
    print("\nCalculating codon weights...")
    weights = relative_adaptiveness(
        sequences=list(reference.values())
    )
    print(
        f"  Calculated weights for {len(weights)} codons."
    )

    # Read GFF
    print("\nReading genome GFF...")
    gff = parse_gff(args.gff)
    print(
        f"  Found {len(gff)} CDS features."
    )

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

    # Calculate statistics
    print(
        "\nCalculating CAI, GC, and GC3s..."
    )
    results = []

    for gene_id, sequence in genome.items():
        if gene_id not in gff:
            print(
                f"WARNING: {gene_id} not found in GFF"
            )
            continue
        info = gff[gene_id]
        try:
            cai = CAI(
                sequence,
                weights=weights
            )
        except Exception as e:
            print(
                f"WARNING: CAI failed for "
                f"{gene_id}: {e}"
            )
            cai = float("nan")
        gc = calculate_gc(sequence)
        gc3s = calculate_gc3s(sequence)
        midpoint = (
            info["start"] + info["end"]
        ) / 2
        results.append({
            "gene_id": gene_id,
            "contig": info["contig"],
            "start": info["start"],
            "end": info["end"],
            "strand": info["strand"],
            "position": midpoint,
            "CAI": cai,
            "GC": gc,
            "GC3s": gc3s,
            "interesting":
                gene_id in interesting_ids
        })


    # Write table
    print(
        "\nWriting results..."
    )
    with open(
        args.output,
        "w",
        newline=""
    ) as outfile:
        writer = csv.writer(
            outfile,
            delimiter="\t"
        )
        writer.writerow([
            "gene_id",
            "contig",
            "start",
            "end",
            "strand",
            "position",
            "coreCAI",
            "GC",
            "GC3s",
            "interesting"
        ])
        for r in results:
            def fmt(x):
                if math.isnan(x):
                    return "NA"
                return f"{x:.6f}"
            writer.writerow([
                r["gene_id"],
                r["contig"],
                r["start"],
                r["end"],
                r["strand"],
                f"{r['position']:.1f}",
                fmt(r["CAI"]),
                fmt(r["GC"]),
                fmt(r["GC3s"]),
                r["interesting"]
            ])
    print(
        f"Saved: {args.output}"
    )


if __name__ == "__main__":
    main()