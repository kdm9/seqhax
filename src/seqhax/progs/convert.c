/*******************************************************************************
*                   convert -- Convert FASTA between FASTQ                    *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
convert_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax convert [options] FILE\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -a     Output FASTA.\n");
    fprintf(stream, "    -q     Output FASTQ (adding qualities).\n");
}
static const char *convert_optstr = "aq";



/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
convert_main(int argc, char *argv[])
{
    int c = 0;
    bool fasta = true;
    const char *filename = NULL;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, convert_optstr)) >= 0) {
        switch (c) {
        case 'a':
            fasta = true;
            break;
        case 'q':
            fasta = false;
            break;
        case 'h':
            convert_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            convert_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (optind >= argc) {
        convert_usage(stderr);
        return EXIT_FAILURE;
    }
    filename = argv[optind];

    struct qes_seq *seq = qes_seq_create();
    struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
    ssize_t res;
    char *tmpqual = NULL;
    size_t tmpqual_sz = 0;
    while ((res = qes_seqfile_read(sf, seq)) > 0) {
        if (fasta) {
            qes_seq_print(seq, stdout, true, 0);
        } else {
            size_t seqlen = seq->seq.len;
            if (seq->qual.len < seqlen) {
                if (seqlen > tmpqual_sz) {
                    tmpqual_sz = seqlen + 1;
                    tmpqual = realloc(tmpqual, seqlen);
                    if (tmpqual == NULL) {
                        fputs("OUT OF MEM\n", stderr);
                        return EXIT_FAILURE;
                    }
                }
                memset(tmpqual, 'I', seqlen);
                tmpqual[seqlen] = 0;
                qes_seq_fill_qual(seq, tmpqual, seqlen);
            }
            qes_seq_print(seq, stdout, false, 0);
        }
    }
    if (tmpqual) free(tmpqual);
    qes_seqfile_destroy(sf);
    qes_seq_destroy(seq);
    return EXIT_SUCCESS;
}
