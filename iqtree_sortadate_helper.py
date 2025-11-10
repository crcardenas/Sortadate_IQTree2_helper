import argparse
import os
import re
import subprocess
from pathlib import Path


def parse_log_file(logfile: str) -> dict:
    """
    Parse the IQ-TREE log file to extract ID -> Name mapping
    for logs with columns:
    Subset  Type  Seqs  Sites  Infor  Invar  Model  Name
    """
    id_to_name = {}
    in_table = False
    header_line = "Subset\tType\tSeqs\tSites\tInfor\tInvar\tModel\tName"

    with open(logfile, "r") as f:
        for line in f:
            stripped = line.strip()
            # Detect start of the table
            if stripped == header_line:
                in_table = True
                continue
            # Detect end of section
            if stripped.startswith("Column meanings"):
                break

            if in_table:
                # Split the line by whitespace or tabs
                parts = line.strip().split()
                if len(parts) >= 8 and parts[0].isdigit():
                    idx = int(parts[0])
                    name = parts[-1]  # last column is always the name
                    name = re.sub(r"^p\d+_", "", name)
                    id_to_name[idx] = name

    return id_to_name

def run_pxrr(treefile: Path, ranked_ogs: str, out_file: Path):
    # assign pxrr command to reroot trees (suppress stdout)
    cmd = [
        "pxrr", "-t", str(treefile), "-r", "-g", ranked_ogs, "-o", str(out_file),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def split_trees(
    treefile: str,
    id_to_name: dict,
    outdir: str,
    suffix: str,
    fasta_dir: str = None,
    reroot: str = None,
    keep_only_outgroup: bool = False,
    name_replace: str = None,
):
    os.makedirs(outdir, exist_ok=True)

    written = 0
    skipped_missing_fasta = 0
    rerooted = 0
    reroot_failed = 0
    no_outgroups = []
    reroot_fail_list = []

    # prepare outgroup list
    outgroup_list = []
    if reroot:
        outgroup_list = [og.strip() for og in reroot.split(",") if og.strip()]

    # prepare name replacement regex if provided
    name_replace_pattern, name_replace_sub = None, None
    if name_replace:
        try:
            name_replace_pattern, name_replace_sub = name_replace.split(",", 1)
        except ValueError:
            print("[ERROR] Invalid --name_replace format. Use 'pattern,replacement'.")
            exit(1)

    with open(treefile, "r") as f:
        for i, line in enumerate(f, start=1):
            if i not in id_to_name:
                continue

            name = id_to_name[i]

            # apply regex name replacement if requested
            if name_replace_pattern and name_replace_sub:
                name = re.sub(name_replace_pattern, name_replace_sub, name)

            # Apply suffix, e.g. _trimalauto
            basename = f"{name}{suffix}" if suffix else name

            # Check fasta file existence if directory is provided
            if fasta_dir:
                fasta_file = os.path.join(fasta_dir, f"{basename}.fasta")
                if not os.path.exists(fasta_file):
                    skipped_missing_fasta += 1
                    continue

            outfile = Path(outdir) / f"{basename}.treefile"
            # write the original treefile (always)
            with open(outfile, "w") as out:
                out.write(line.strip() + "\n")
            written += 1

            # Handle rerooting if requested
            if reroot:
                if keep_only_outgroup:
                    found = False
                    for og in outgroup_list:
                        # word-boundary match to avoid partial matches
                        if re.search(rf"\b{re.escape(og)}\b", line):
                            found = True
                            break
                    if not found:
                        no_outgroups.append(str(outfile))
                        continue

                rr_file = outfile.with_name(outfile.name + ".rr")

                try:
                    run_pxrr(outfile, reroot, rr_file)
                    rerooted += 1
                except subprocess.CalledProcessError as e:
                    print(f"[ERROR] pxrr failed for {outfile}: {e}")
                    reroot_failed += 1
                    reroot_fail_list.append(str(outfile))

    # Print log summary
    print(f"Trees written: {written}")
    if fasta_dir:
        print(f"Trees skipped (missing fasta): {skipped_missing_fasta}")
    if reroot:
        print(f"Trees rerooted: {rerooted}")
        print(f"Trees reroot failed: {reroot_failed}")
        print(f"Trees without specified outgroups: {len(no_outgroups)}")

    # Always write the no_outgroups list (may be empty per request)
    if reroot:
        list_file = Path(outdir) / "no_outgroups.list"
        with open(list_file, "w") as lf:
            for item in no_outgroups:
                lf.write(item + "\n")
        print(f"Wrote no-outgroups list: {list_file}")

        # If any pxrr failures, write a file listing them
        if reroot_fail_list:
            fail_file = Path(outdir) / "reroot_failed.list"
            with open(fail_file, "w") as ff:
                for item in reroot_fail_list:
                    ff.write(item + "\n")
            print(f"Wrote reroot-failures list: {fail_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Split multi-tree IQ-TREE output into per-locus tree files, with optional rerooting."
    )
    parser.add_argument(
        "-l", "--logfile", required=True,
        help="IQ-TREE log file with locus table; required",
    )
    parser.add_argument(
        "-i", "--input_locus_trees", required=True,
        help="Multi-tree .treefile from IQ-TREE; required",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory for individual trees; required",
    )
    parser.add_argument(
        "-s",  "--suffix",  required=False,
        help="Optional suffix for file names (default: '_trimalauto')",
    )
    parser.add_argument(
        "-f",  "--fasta_directory", required=False,
        help="Optional directory of fasta files to check existence before writing tree",
    )
    parser.add_argument(
        "-r", "--reroot", type=str, required=False,
        help="Comma-separated outgroup list for rerooting (will be passed to pxrr -g)",
    )
    parser.add_argument(
        "--keep_only_outgroup", action="store_true", required=False,
        help="If set, only reroot trees containing at least one specified outgroup. Others are listed in no_outgroups.list; requires -r/--reroot flag",
    )
    parser.add_argument(
        "--name_replace", type=str, required=False,
        help=r"Optional regex replacement applied to locus names before suffix (format: 'pattern,replacement'): e.g., \"uce(\\d+),uce-\\1.\".",
    )

    args = parser.parse_args()

    # assign function logic
    id_to_name = parse_log_file(args.logfile)
    split_trees(
        args.input_locus_trees,
        id_to_name,
        args.output,
        args.suffix,
        args.fasta_directory,
        args.reroot,
        args.keep_only_outgroup,
        args.name_replace,
    )


if __name__ == "__main__":
    main()
