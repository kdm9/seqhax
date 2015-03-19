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
            print >> sys.stderr, "WARNING: truncating", seq.name
            seq_len = qual_len
        print "@{}\n{}\n+\n{}".format(seq.name,
                                      seq.sequence[:seq_len],
                                      qual.accuracy[:seq_len])

if __name__ == "__main__":
    main()
