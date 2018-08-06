/*******************************************************************************
*                   stats -- Basic stats on sequence files                    *
*******************************************************************************/


#include <omp.h>
#include <iostream>
#include <string>
#include <cctype>
#include <stdexcept>
#include <atomic>



#include "gzstream.h"

using namespace std;

static void
stats_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax stats FILE ...\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -t THREADS  Number of parallel jobs [1]\n");
    fprintf(stream, "\nFILEs must be FASTQs, optionally gzip-compressed\n");
}
static const char *stats_optstr = "t:";



/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
stats_main(int argc, char *argv[])
{
    int c = 0;
    int n_threads = 1;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, stats_optstr)) >= 0) {
        switch (c) {
        case 't':
            n_threads = atoi(optarg);
            break;
        case 'h':
            stats_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            stats_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (optind >= argc) {
        stats_usage(stderr);
        return EXIT_FAILURE;
    }

    atomic_size_t errorcount(0);
    cout << "filename\treads\tacgtnbases\tacgtbases\tgcbases\tstatus" << endl;
    #pragma omp parallel for num_threads(n_threads) schedule(dynamic) shared(errorcount)
    for (int file = optind; file < argc; file++) {
        size_t nbases = 0;
        size_t ngc = 0;
        size_t nnon = 0;
        size_t nlines = 0;
        size_t seqlen = 0;
        string status = "OK";

        igzstream fp;
        try {
            fp.open(argv[file]);
            string line;
            for (; getline(fp, line); nlines++) {
                switch (nlines % 4) {
                    case 0:
                        if (line[0] != '@') throw runtime_error("Header doesn't start with @");
                        break;
                    case 1:
                        for (size_t i = 0; i < line.size(); i++) {
                            switch (tolower(line[i])) {
                                case 'c':
                                case 'g':
                                    ngc++;
                                case 'a':
                                case 't':
                                    nnon++;
                                case 'n':
                                    nbases++;
                                    break;
                                default:
                                    throw runtime_error("Non-ACGTN in sequence: " + line[i]);
                            }
                        }
                        seqlen = line.size();
                        break;
                    case 2:
                        if (line[0] != '+') throw runtime_error("Qualheader doesn't start with +");
                        break;
                    case 3:
                        if (line.size() != seqlen) throw runtime_error("Size of sequence (" + to_string(seqlen) + ") and quality (" + to_string(line.size()) +") doesn't match");
                        break;
                }
            }
            if (nlines == 0) throw runtime_error("Empty file");
            if (nlines % 4 != 0) throw runtime_error("Extra lines in file (not groups of four)");
        } catch (exception &e) {
            status = e.what();
            errorcount++;
        }

        #pragma omp critical
        {
            cout << argv[file] << '\t' << (nlines/4) << '\t' << nbases << '\t' << nnon << '\t' << ngc << '\t' << status << endl;
        }
    }
    return errorcount == 0 ? 0 : 1;
}

