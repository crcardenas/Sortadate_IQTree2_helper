# Cody Raul Cardenas
# September 27th 2025
# Python 3.10
# requires phyx 1.1 in your environment
import argparse
import os
import re
import subprocess
from pathlib import Path


def parse_log_file(logfile: str) -> dict:
# Parse the IQ-TREE log file to extract ID -> Name mapping.
# Stops at the 'Column meanings' section.
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


def run_pxrr(treefile: Path, ranked_ogs: str, out_file: Path):
# Run pxrr to reroot the trees.
# Parameters
#    treefile : Path
#        Input tree file to reroot
#    ranked_ogs : str
#        Comma-separated outgroup list (or whatever pxrr expects for -g)
#    out_file : Path
#        Path to write the rerooted tree (pxrr -o)
    cmd = ["pxrr", "-t", str(treefile), "-r", "-g", ranked_ogs, "-o", str(out_file),]
    print(f"[RUN] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def split_trees(
    treefile: str,
    id_to_name: dict,
    outdir: str,
    suffix: str,
    fasta_dir: str = None,
    reroot: str = None,
    keep_only_outgroup: bool = False,
):
    os.makedirs(outdir, exist_ok=True)

    written = 0
    skipped_missing_fasta = 0
    rerooted = 0
    reroot_failed = 0
    no_outgroups = []
    reroot_fail_list = []

    # prepare outgroup list (trim whitespace)
    outgroup_list = []
    if reroot:
        outgroup_list = [og.strip() for og in reroot.split(",") if og.strip()]

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
                    skipped_missing_fasta += 1
                    continue

            outfile = Path(outdir) / f"{basename}.treefile"
            # write the original treefile (always)
            with open(outfile, "w") as out:
                out.write(line.strip() + "\n")
            written += 1

            # Handle rerooting if requested
            if reroot:
                # If keep_only_outgroup is set, check that at least one outgroup is present in the tree string
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

                # Create rerooted filename by appending .rr to the original filename
                # e.g. core_10000_trimalauto.treefile -> core_10000_trimalauto.treefile.rr
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
        "-l", "--logfile", required=True, help="IQ-TREE log file with locus table",
    )
    parser.add_argument( "-i","--input_locus_trees",required=True,help="Multi-tree .treefile from IQ-TREE",
    )
    parser.add_argument(
        "-o","--output",required=True,help="Output directory for individual trees",
    )
    parser.add_argument(
        "-s","--suffix",required=False,default="_trimalauto",help="Optional suffix for file names (default: '_trimalauto')",
    )
    parser.add_argument(
        "-f","--fasta_directory",required=False,help="Optional directory of fasta files to check existence before writing tree",
    )
    parser.add_argument(
        "-r","--reroot",type=str,required=False,help="Comma-separated outgroup list for rerooting (will be passed to pxrr -g)",
    )
    parser.add_argument(
        "--keep_only_outgroup",action="store_true", help="If set, only reroot trees containing at least one specified outgroup. Others are listed in no_outgroups.list",
    )

    args = parser.parse_args()

    id_to_name = parse_log_file(args.logfile)
    split_trees(
        args.input_locus_trees,
        id_to_name,
        args.output,
        args.suffix,
        args.fasta_directory,
        args.reroot,
        args.keep_only_outgroup,
    )


if __name__ == "__main__":
    main()

    main()
