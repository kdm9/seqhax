#include "kmseq.hh"
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
        size_t count = 0;
        for (KSeq s; seqs.next_read(s); count++) {
            cout << s;
        }
    }
    return 0;
}
