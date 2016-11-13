/*******************************************************************************
*                   stats -- Basic stats on sequence files                    *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
stats_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax stats FILE ...\n");
}
static const char *stats_optstr = "";



/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
stats_main(int argc, char *argv[])
{
    int c = 0;
    const char *filename = NULL;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, stats_optstr)) >= 0) {
        switch (c) {
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

    puts("filename\treads\tbases");
    for (; optind < argc; optind++) {
        filename = argv[optind];
        struct qes_seq *seq = qes_seq_create();
        struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
        size_t n_reads = 0;
        size_t n_bp = 0;
        while (qes_seqfile_read(sf, seq) > 0) {
            n_reads++;
            n_bp += seq->seq.len;
        }
        if (n_reads == 0) {
            fprintf(stderr, "WARNING: Invalid or empty file '%s'\n", filename);
        }
        printf("%s\t%zu\t%zu\n", filename, n_reads, n_bp);
        qes_seqfile_destroy(sf);
        qes_seq_destroy(seq);
    }
    return EXIT_SUCCESS;
}
