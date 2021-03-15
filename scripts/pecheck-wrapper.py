#!/usr/bin/env python3
import argparse as ap
import subprocess as sp
from collections import defaultdict
import fnmatch
from os import path as op
import os
import re
import multiprocessing as mp
from sys import stdout, stderr


def run(cmd):
    return (cmd, sp.check_call(cmd, shell=True, executable="/bin/bash"))


def main(inputs, outdir, dryrun=True, jobs=1, force=False, gzip="gzip", fastqpattern='*.fastq.gz', rgpattern="_R([12])", lane_pattern=r"_L\d+", expect_n=2):
    if not op.exists(outdir) and not dryrun:
        os.makedirs(outdir)

    fastqs = []
    for fileordir in inputs:
        if op.isfile(fileordir) and fnmatch.fnmatch(fileordir, fastqpattern):
            fastqs.append(fileordir)
        else:
            for root, dirs, files in os.walk(fileordir):
                for file in fnmatch.filter(files, fastqpattern):
                    fastqs.append(op.join(root, file))

    matching_fastqs = defaultdict(list)
    for fastq in fastqs:
        group = re.sub(rgpattern, "_READGROUPPATTERN", fastq)
        group = re.sub(lane_pattern, "_LANEPATTERN", group)
        matching_fastqs[group].append(fastq)

    todo = []
    for group, members in matching_fastqs.items():
        if len(members) != expect_n:
            print("ERROR: group of matching files does not contain", expect_n, "files:", group, members, file=stderr)
            continue
        out = op.basename(group).replace("_READGROUPPATTERN", "_il").replace("_LANEPATTERN", "")
        if out.endswith(".fq") or out.endswith(".fq.gz"):
            out = re.sub(r"\.fq(\.gz)?$", ".fastq.gz", out)
        if not out.endswith(".gz"):
            out += ".gz"
        out = op.join(outdir, out)
        tmpout = out + "_TMPOUT"
        fqs = "' '".join(sorted(members))
        if force:
            ifnotexist = ""
        else:
            ifnotexist = f"test -f '{out}' || "
        cmd = f"({ifnotexist} ( seqhax pecheck -o >( {gzip} > '{tmpout}' ) '{fqs}' >'{out}.tsv' && mv '{tmpout}' '{out}' )) 2>{out}.log"
        todo.append(cmd)

    if dryrun:
        if len(todo) > 0:
            print(*todo, sep="\n")
        else:
            print("No groups found", file=stderr)
        return
    pool = mp.Pool(jobs)
    for cmd, res in pool.imap(run, todo):
        print(cmd)
        if isinstance(res, sp.CalledProcessError):
            pool.terminate()
            raise res


if __name__ == "__main__":
    a = ap.ArgumentParser("pecheck-wrapper.py")
    a.add_argument("--gzip", type=str, default="gzip",
                   help="Which command should be used to gzip? (try pigz!)")
    a.add_argument("-L", "--lane-pattern", default=r"_L\d+(?=[_.])", type=str,
                   help="regex to match the lane number (typcally L1 or L001)")
    a.add_argument("-R", "--readgroup-pattern", default=r"_R([12](?=[_.])", type=str,
                   help="regex to match the read group (typcally R1/R2, sometimes just _1/_2)")
    a.add_argument("-F", "--fastq-pattern", default="*.fastq.gz", type=str,
                   help="shell glob pattern to match a fastq file.")
    a.add_argument("-j", "--jobs", default=1, type=int,
                   help="number of parallel jobs")
    a.add_argument("-f", "--force", action="store_true",
                   help="Force creation of merged outputs even if they exist")
    a.add_argument("-N", "--expect-n", default=2, type=int,
                   help="Expect N files per group.")
    a.add_argument("-n", "--dry-run", action="store_true",
                   help="Don't actually do anything, just print what would be done")
    a.add_argument("OUTPUT", type=str,
                   help="Output directory (created if non-existant)")
    a.add_argument("INPUT", nargs='+', type=str,
                   help="Input directory(s) or files")
    args = a.parse_args()
    main(args.INPUT, args.OUTPUT, dryrun=args.dry_run, jobs=args.jobs,
         force=args.force, gzip=args.gzip, fastqpattern=args.fastq_pattern,
         rgpattern=args.readgroup_pattern, lane_pattern=args.lane_pattern, expect_n=args.expect_n)
