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
    print(cmd)
    return sp.check_call(cmd, shell=True, executable="/bin/bash"),


def main(inputs, outdir, dryrun=True, jobs=1, force=False, gzip="gzip", fastqpattern='*.fastq.gz'):
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
        sans_r12 = re.sub("_R[12]", "", fastq)
        matching_fastqs[sans_r12].append(fastq)

    todo = []
    for group, members in matching_fastqs.items():
        if len(members) != 2:
            print("ERROR: group of matching files does not contain a pair:", group, members, file=stderr)
            continue
        out = op.basename(group)
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
    res = pool.imap(run, todo)
    for r in res:
        if isinstance(r, sp.CalledProcessError):
            pool.terminate()
            raise r


if __name__ == "__main__":
    a = ap.ArgumentParser("pecheck-wrapper.py")
    a.add_argument("--gzip", type=str, default="gzip",
                   help="Which command should be used to gzip? (try pigz!)")
    a.add_argument("-j", "--jobs", default=1, type=int,
                   help="number of parallel jobs")
    a.add_argument("-f", "--force", action="store_true",
                   help="Force creation of merged outputs even if they exist")
    a.add_argument("-n", "--dry-run", action="store_true",
                   help="Don't actually do anything, just print what would be done")
    a.add_argument("OUTPUT", type=str,
                   help="Output directory (created if non-existant)")
    a.add_argument("INPUT", nargs='+', type=str,
                   help="Input directory(s) or files")
    args = a.parse_args()
    main(args.INPUT, args.OUTPUT, dryrun=args.dry_run, jobs=args.jobs,
         force=args.force, gzip=args.gzip)
