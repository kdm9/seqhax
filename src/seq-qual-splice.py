#!/usr/bin/env python
from __future__ import print_function
import screed
import sys
from itertools import izip

def main():
    seq_file = sys.argv[1]
    qual_file = sys.argv[2]

    seqs = screed.open(seq_file)
    quals = screed.open(qual_file)

    for seq, qual in izip(seqs, quals):
        seq_len = len(seq.sequence)
        qual_len = len(qual.accuracy)
        if qual_len < seq_len:
            print("WARNING: truncating", seq.name, file=sys.stderr)
            seq_len = qual_len
        out_seq = "@{}\n{}\n+\n{}\n".format(seq.name, seq.sequence[:seq_len],
                                          qual.accuracy[:seq_len])
        print(out_seq, end='')

USAGE = """
USAGE: {} <seq_fastq> <qual_fastq>

Splices the sequence lines from <seq_fastq> with the quality lines of
<qual_fastq>, truncating both the the shorter of the two.
"""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(USAGE.format(sys.argv[0]), file=sys.stderr)
        exit()
    main()
