# Cody Raul Cardenas
# September 27th 2025
# Python 3.10

import argparse
import os
import re


def parse_log_file(logfile):
#    Parse the IQ-TREE log file to extract ID -> Name mapping.
#    Stops at the 'Column meanings' section.
    id_to_name = {}
    with open(logfile, "r") as f:
        for line in f:
            if line.strip().startswith("Column meanings"):
                break
            match = re.match(
                r"\s*(\d+)\s+DNA\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\S+\s+(\S+)",
                line,
            )
            if match:
                idx, name = match.groups()
                # Strip off prefix like p123_
                name = re.sub(r"^p\d+_", "", name)
                id_to_name[int(idx)] = name
    return id_to_name


def split_trees(treefile, id_to_name, outdir, suffix, fasta_dir=None):
    os.makedirs(outdir, exist_ok=True)

    written = 0
    skipped = 0

    with open(treefile, "r") as f:
        for i, line in enumerate(f, start=1):
            if i not in id_to_name:
                continue

            name = id_to_name[i]

            # Apply suffix, e.g. _trimalauto
            basename = f"{name}{suffix}" if suffix else name

            # Check fasta file existence if directory is provided
            if fasta_dir:
                fasta_file = os.path.join(fasta_dir, f"{basename}.fasta")
                if not os.path.exists(fasta_file):
                    skipped += 1
                    continue

            outfile = os.path.join(outdir, f"{basename}.treefile")
            with open(outfile, "w") as out:
                out.write(line.strip() + "\n")
            written += 1

    # Print log summary
    print(f"Trees written: {written}")
    if fasta_dir:
        print(f"Trees skipped (missing fasta): {skipped}")


def main():
    parser = argparse.ArgumentParser(
        description="Split multi-tree IQ-TREE output into per-locus tree files."
    )
    parser.add_argument(
        "-l", "--logfile", required=True, help="IQ-TREE log file with locus table"
    )
    parser.add_argument(
        "-i", "--input_locus_trees", required=True, help="Multi-tree .treefile from IQ-TREE"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output directory for individual trees"
    )
    parser.add_argument(
        "-s", "--suffix", required=False, help="Optional suffix for file names (default: '_trimalauto')",
    )
    parser.add_argument(
        "-f", "--fasta_directory", required=False, help="Optional directory of fasta files to check existence before writing tree",
    )

    args = parser.parse_args()

    id_to_name = parse_log_file(args.logfile)
    split_trees(
        args.input_locus_trees,
        id_to_name,
        args.output,
        args.suffix,
        args.fasta_directory,
    )


if __name__ == "__main__":
    main()
