#include "kmseq.hh"
#include <vector>
#include <iostream>
using namespace std;
using namespace kmseq;


int cat_main(int argc, char *argv[])
{
    if (argc < 2) {
        cerr << "USAGE: seqhax cat FILE...\n";
        return 1;
    }
    for (int i = 1; i < argc; i++) {
        KSeqReader seqs(argv[i]);
        vector<KSeq> sv;
        while (seqs.next_chunk(sv, 1000) > 0) {
            for (const auto &s: sv) {
                cout << s;
            }
        }
    }
    return 0;
}
