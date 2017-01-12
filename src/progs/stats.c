/*******************************************************************************
*                   stats -- Basic stats on sequence files                    *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>

#include "qes_seqfile.h"


static void
stats_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax stats FILE ...\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -t THREADS  Number of parallel jobs [1]\n");
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

    puts("filename\treads\tbases\tperc_gc\tnon_acgt");
    #pragma omp parallel for num_threads(n_threads) schedule(dynamic)
    for (int file = optind; file < argc; file++) {
        const char *filename = argv[file];
        struct qes_seq *seq = qes_seq_create();
        struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
        size_t n_reads = 0;
        size_t n_bp = 0;
        size_t comp[128] = {0};
        while (qes_seqfile_read(sf, seq) > 0) {
            n_reads++;
            n_bp += seq->seq.len;
            for (size_t i = 0; i < seq->seq.len; i++) {
                comp[seq->seq.str[i] & 0xdf] += 1;
            }
        }
        ssize_t at = comp['A'] + comp['T'];
        ssize_t gc = comp['C'] + comp['G'];
        ssize_t acgt = at + gc;
        ssize_t non_acgt = n_bp - acgt;
        float prop_gc = (float) gc / (float) acgt;
        #pragma omp critical
        {
            if (n_reads == 0) {
                fprintf(stderr, "WARNING: Invalid or empty file '%s'\n", filename);
            }
            printf("%s\t%zu\t%zu\t%0.1f\t%zu\n", filename, n_reads, n_bp,
                        prop_gc * 100, non_acgt);
            fflush(stdout);
        }
        qes_seqfile_destroy(sf);
        qes_seq_destroy(seq);
    }
    return EXIT_SUCCESS;
}
