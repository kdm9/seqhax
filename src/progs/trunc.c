/*******************************************************************************
*                         trunc -- Truncate sequences                         *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
trunc_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax seq [options] FILE\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -l LEN     Fixed (integer) length to truncate at.\n");
    fprintf(stream, "    -p PROP    Truncate at PROP * length bases\n");
    fprintf(stream, "\n");
    fprintf(stream, "At least one of `-l` and `-p` MUST be given.\n");
    fprintf(stream, "FILE should be a sequence file in FASTA or FASTQ format.\n");
    fprintf(stream, "To accept reads from standard input, use '/dev/stdin' as\n");
    fprintf(stream, "the input file.\n");
}
static const char *trunc_optstr = "l:p:";



/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
trunc_main(int argc, char *argv[])
{
    int c = 0;
    size_t truncat = 0;
    float truncprop = 0.;
    const char *filename = NULL;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, trunc_optstr)) >= 0) {
        switch (c) {
        case 'p':
            truncprop = atof(optarg);
            break;
        case 'l':
            truncat = atoi(optarg);
            break;
        case 'h':
            trunc_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            trunc_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (optind >= argc || (truncat == 0 && truncprop == 0.0)) {
        trunc_usage(stderr);
        return EXIT_FAILURE;
    }
    filename = argv[optind];

    struct qes_seq *seq = qes_seq_create();
    struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
    ssize_t res;
    while ((res = qes_seqfile_read(sf, seq)) > 0) {
        if (truncprop != 0.0) {
            truncat = truncprop * seq->seq.len;
        }
        qes_seq_truncate(seq, truncat);
        qes_seq_print(seq, stdout, !qes_seq_has_qual(seq));
        
    }

    qes_seqfile_destroy(sf);
    qes_seq_destroy(seq);
    return EXIT_SUCCESS;
}
