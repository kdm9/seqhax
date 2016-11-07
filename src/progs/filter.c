/*******************************************************************************
*                         filter -- Filters sequences                         *
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
filter_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax filter [options] FILE\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -l LENGTH  Minimum length of each read. [default 1]\n");
    fprintf(stream, "    -f         Output as fasta (no qualities)\n");
    fprintf(stream, "    -p         Paired mode: reads are kept/discared in pairs\n");
    fprintf(stream, "\n");
    fprintf(stream, "FILE should be a sequence file in FASTA or FASTQ format.\n");
    fprintf(stream, "To accept reads from standard input, use '/dev/stdin' as\n");
    fprintf(stream, "the input file.\n");
}
static const char *filter_optstr = "pfhl:";

/*******************************************************************************
*                                   FILTERS                                   *
*******************************************************************************/

typedef struct {
    bool r1;
    bool r2;
}  filter_pass_t;


void
length_filter(filter_pass_t *pass, struct qes_seq *seq1, struct qes_seq *seq2,
              int min_len)
{
    if (pass->r1 && seq1->seq.len < min_len) {
        pass->r1 = false;
    }
    if (pass->r2 && seq2->seq.len < min_len) {
        pass->r2 = false;
    }
}


/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
filter_main(int argc, char *argv[])
{
    int c = 0;
    int min_len = 0;
    bool fasta = false;
    bool paired = false;
    const char *filename = NULL;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, filter_optstr)) >= 0) {
        switch (c) {
        case 'f':
            fasta = true;
            break;
        case 'p':
            paired = true;
            break;
        case 'l':
            min_len = (size_t) atol(optarg);
            break;
        case 'h':
            filter_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            filter_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (optind >= argc) {
        filter_usage(stderr);
        return EXIT_FAILURE;
    }
    filename = argv[optind];

    struct qes_seq *r1 = qes_seq_create();
    struct qes_seq *r2 = qes_seq_create();
    struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
    ssize_t res1, res2;
    while ((res1 = qes_seqfile_read(sf, r1)) > 0) {
        filter_pass_t pass = {true, true};
        res2 = qes_seqfile_read(sf, r2);
        if (res2 <= 0) {
            // In case we have an odd number of reads
            pass.r2 = false;
        }

        if (min_len > 0) {
            length_filter(&pass, r1, r2, min_len);
        }

        if (paired) {
            if (pass.r1 && pass.r2) {
                qes_seq_print(r1, stdout, fasta, 1);
                qes_seq_print(r2, stdout, fasta, 2);
            }
        } else {
            if (pass.r1) {
                qes_seq_print(r1, stdout, fasta, 0);
            }
            if (pass.r2) {
                qes_seq_print(r2, stdout, fasta, 0);
            }
        }
    }

    qes_seqfile_destroy(sf);
    qes_seq_destroy(r1);
    qes_seq_destroy(r2);
    return EXIT_SUCCESS;
}
