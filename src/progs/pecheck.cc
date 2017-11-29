#include <iostream>
#include <fstream>
#include <iterator>
#include <string>

#include <getopt.h>

#include "kmseq.hh"
#include "libseqhax.hh"

using namespace std;
using namespace kmseq;
using namespace seqhax;

struct PECheckOptions {
    string r1file;
    string r2file;
    string outfile;
};

inline void helpmsg(void)
{
    cerr << "USAGE:" << endl;
    cerr << "    seqhax pecheck [-o OUTPUT] r1 r2" << endl
         << endl;
    cerr << "OPTIONS:"<< endl;
    cerr << "    -o FILE    Output interleaved reads to FILE. Use - for stdout. (default: no output)"<< endl;
}

int parse_args(PECheckOptions &opt, int argc, char *argv[])
{

    int c;
    int ret = EXIT_SUCCESS;
    while ((c = getopt(argc, argv, "o:")) > 0) {
        switch (c) {
            case 'o':
                opt.outfile = optarg;
                if (opt.outfile == "-") opt.outfile = "/dev/fd/1";
                break;
            case '?':
                ret = EXIT_FAILURE;
            case 'h':
                helpmsg();
                return ret;
        }
    }

    if (optind + 1 >= argc) {
        helpmsg();
        return EXIT_FAILURE;
    }

    opt.r1file = argv[optind++];
    opt.r2file = argv[optind++];
    return 0;
}


int pecheck_main(int argc, char *argv[])
{
    PECheckOptions opt;
    ofstream output;

    if (parse_args(opt, argc, argv) != 0) {
        return 1;
    }
    cerr << "From " << opt.r1file << " and " << opt.r2file << " to " << opt.outfile << endl;

    if (opt.outfile != "") {
        output.open(opt.outfile);
    }

    KSeqPairReader seqs(opt.r1file, opt.r2file);
    size_t count = 0;
    for (KSeqPair sp; seqs.next_pair(sp); count++) {
        string n1, n2;
        extract_readname(sp.r1.name, n1);
        extract_readname(sp.r2.name, n2);
        if (n1 != n2) {
            cerr << "Pair mismatch! Are you sure these are paired-end reads?" << endl;
            cerr << endl;
            cerr << "At pair " << count << endl;
            cerr << "'" << n1 << "' != '" << n2 << "'" << endl;
            return 1;
        }
        if (output.is_open()) {
            output << sp;
        }
    }
    cerr << "Processed " << count << " pairs" << endl;
    return 0;
}

#ifdef PECHECK_STANDALONE
int main(int argc, char *argv[])
{
    return pecheck_main(argc, argv);
}
#endif
