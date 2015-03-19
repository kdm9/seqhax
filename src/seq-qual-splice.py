#!/usr/bin/env python
# Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
