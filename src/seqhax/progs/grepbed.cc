#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <regex>
#include "kmseq.hh"


using namespace std;
using namespace kmseq;


inline void grepbed_helpmsg(void)
{
    cerr << "USAGE:" << endl;
    cerr << "    seqhax grepbed -p REGEX FASTA_FILE" << endl
         << endl;
    cerr << "OPTIONS:"<< endl;
    cerr << "    -p REGEX   A nucleotide regex to search for. E.g. CHH = C[ACT][ACT]."<< endl;
}


int grepbed_main(int argc, char *argv[])
{
    std::string pattern;
    std::string filename;
    int ret = EXIT_SUCCESS;

    int c;
    while ((c = getopt(argc, argv, "p:h")) > 0) {
        switch (c) {
            case 'p':
                pattern = optarg;
                break;
            case '?':
                ret = EXIT_FAILURE;
            case 'h':
                grepbed_helpmsg();
                return ret;
        }
    }
    if (optind >= argc) {
        grepbed_helpmsg();
        return EXIT_FAILURE;
    }

    filename = argv[optind++];

    KSeqReader seqs(filename);
    
    for (KSeq seq; seqs.next_read(seq);) {
        regex re(string("(?=(") + pattern + string("))."));
        sregex_iterator next(seq.seq.begin(), seq.seq.end(), re);
        sregex_iterator end;
        for (; next != end; *next++) {
            smatch match = *next;
            cout << seq.id() << "\t" << match.position() + 1 << "\t" << match.str(1) << endl;
        }
    }
    return ret;
}

#ifdef SEQHAX_STANDALONE
int main(int argc, char *argv[])
{
    return grepbed_main(argc, argv);
}
#endif

