#include <iostream>
#include <fstream>
#include <iterator>
#include <string>

#include <getopt.h>
#include <omp.h>

#include "kmseq.hh"
#include "libseqhax.hh"

using namespace std;
using namespace kmseq;
using namespace seqhax;

struct PECheckOptions {
    vector<string> inputfiles;
    string outfile;
    bool print_table;
    int num_threads;
    bool interleaved;
};

inline void helpmsg(void)
{
    cerr << "USAGE:" << endl
         << "    seqhax pecheck [OPTIONS] r1fq r2fq [r1fq r2fq ...]" << endl
         << "    seqhax pecheck [OPTIONS] -i ilfq [ilfq ilfq ...]" << endl
         << endl
         << "OPTIONS:"<< endl
         << "    -o FILE        Output interleaved reads to FILE. Use - for stdout." << endl
         << "                   All sets of paired end files will end up in same" << endl
         << "                   output!!! (default: no output)" << endl
         << "    -i             Interleaved inputs" << endl
         << "    -q             Don't print tabular summary" << endl
         << "    -t THREADS     Number of parallel threads (default: no output)"<< endl;
}

int parse_args(PECheckOptions &opt, int argc, char *argv[])
{

    int c;
    int ret = EXIT_SUCCESS;
    opt.print_table = true;
    opt.num_threads = 1;
    opt.interleaved = false;
    while ((c = getopt(argc, argv, "o:t:i")) > 0) {
        switch (c) {
            case 'o':
                opt.outfile = optarg;
                if (opt.outfile == "-") {
                    opt.outfile = "/dev/fd/1";
                    opt.print_table = false;
                }
                break;
            case 't':
                opt.num_threads = atoi(optarg);
                break;
            case 'i':
                opt.interleaved = true;
                break;
            case 'q':
                opt.print_table = false;
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
    for (; optind + 1 < argc;) {
        opt.inputfiles.push_back(argv[optind++]);
    }
    if (opt.inputfiles.size() > 1 && opt.outfile != "") {
        cerr << "WARNING: combining multiple pairs of files to a single output. This might be a **very** bad idea." << endl;
        opt.num_threads = 1;
    }
    return 0;
}


int pecheck_main(int argc, char *argv[])
{
    PECheckOptions opt;
    ofstream output;

    if (parse_args(opt, argc, argv) != 0) {
        return 1;
    }

    if (opt.outfile != "") {
        cerr << "Interleaving reads to " << opt.outfile << endl;
        output.open(opt.outfile);
        opt.num_threads = 1;
    }

    if (opt.print_table) {
        cout << "r1_file\tr2_file\tstatus\tread_pairs\n";
    }
    bool allpass = true;
    #pragma omp parallel for schedule(dynamic) num_threads(opt.num_threads) shared(allpass)
    for (size_t i = 0; i < opt.inputfiles.size(); i++) {
        KSeqPairReader seqs;
        try {
            if (opt.interleaved) {
                seqs.open(opt.inputfiles[i]);
            } else {
                seqs.open(opt.inputfiles[i], opt.inputfiles[i+1]); i++;
            }
        } catch (exception &exc) {
            cerr << exc.what() << endl;
            break;
        }
                
        size_t count = 0;
        bool pass = true;
        for (KSeqPair sp; seqs.next_pair(sp); count++) {
            string n1, n2;
            int mate1, mate2;
            extract_readname(sp.r1.name, n1, mate1);
            extract_readname(sp.r2.name, n2, mate2);
            if (n1 != n2 || mate1 != 1 || mate2 != 2) {
                allpass = pass = false;
                #pragma omp critical
                {
                    cerr << "Pair mismatch! Are you sure these are paired-end reads?" << endl;
                    cerr << endl;
                    cerr << "At pair " << count << endl;
                    cerr << "'" << n1 << "' != '" << n2 << "'" << endl;
                    cerr << mate1 << " != " << mate2 << endl;
                }
                break;
            }
            if (output.is_open()) {
                // This won't need a critical section, as we disable multi-threading when printing reads.
                output << sp;
            }
        }
        if (opt.print_table) {
            #pragma omp critical
            {
                string status = pass ? "OK" : "FAIL";
                cout << opt.inputfiles[i-2] << "\t" << opt.inputfiles[i-1] << "\t"
                     << status << "\t" << count << endl;
            }
        }
    }
    if (allpass) {
        cerr << "All sets of reads were correctly paired" << endl;
        return EXIT_SUCCESS;
    } else {
        cerr << "Not all reads were correctly paired" << endl;
        return EXIT_FAILURE;
    }
}

#ifdef PECHECK_STANDALONE
int main(int argc, char *argv[])
{
    return pecheck_main(argc, argv);
}
#endif
