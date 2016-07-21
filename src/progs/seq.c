/*******************************************************************************
*                      seq -- Misc sequence modification                      *
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
seq_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax seq [options] FILE\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -P SEQ     Add SEQ as prefix, adding quality charachters if fastq.\n");
    fprintf(stream, "    -S SEQ     Add SEQ as suffix, adding quality charachters if fastq.\n");
    fprintf(stream, "    -p         Paired mode: reads are kept/discared in pairs\n");
    fprintf(stream, "\n");
    fprintf(stream, "FILE should be a sequence file in FASTA or FASTQ format.\n");
    fprintf(stream, "To accept reads from standard input, use '/dev/stdin' as\n");
    fprintf(stream, "the input file.\n");
}
static const char *seq_optstr = "fP:S:";


int
seq_print_presuf(struct qes_seq *seq, FILE *stream,
                 char *prefix, size_t prefixlen,
                 char *suffix, size_t suffixlen,
                 bool fasta)
{
    size_t iii = 0;
    if (stream == NULL || seq == NULL) {
        return 1;
    }

    if (seq->qual.len < 1) {
        fasta = true;
    }

    if (fasta) {
        fputc('>', stream);
    } else {
        fputc('@', stream);
    }
    fputs(seq->name.str, stream);
    if (seq->comment.str) {
        fputc(' ', stream);
        fputs(seq->comment.str, stream);
    }
    fputc('\n', stream);
    if (prefix != NULL) {
        fputs(prefix, stream);
    }
    fputs(seq->seq.str, stream);
    if (suffix != NULL) {
        fputs(suffix, stream);
    }
    fputc('\n', stream);
    if (!fasta) {
        fputs("+\n", stream);
        /* qual score of I is valid for all encodings */
        for (iii = 0; prefix != NULL &&  iii < prefixlen; iii ++) {
            fputc('I', stream);
        }
        fputs(seq->qual.str, stream);
        for (iii = 0; suffix != NULL && iii < suffixlen; iii ++) {
            fputc('I', stream);
        }
        fputc('\n', stream);
    }
    return 0;
}

/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
seq_main(int argc, char *argv[])
{
    int c = 0;
    bool fasta = false;
    char *prefix = NULL;
    size_t prefixlen = 0;
    char *suffix = NULL;
    size_t suffixlen = 0;
    const char *filename = NULL;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, seq_optstr)) >= 0) {
        switch (c) {
        case 'f':
            fasta = true;
            break;
        case 'P':
            prefix = strdup(optarg);
            prefixlen = strlen(prefix);
            break;
        case 'S':
            suffix = strdup(optarg);
            suffixlen = strlen(suffix);
            break;
        case 'h':
            seq_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            seq_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (fasta) {
        fprintf(stderr, "fasta\n");
    }
    if (optind >= argc) {
        seq_usage(stderr);
        return EXIT_FAILURE;
    }
    filename = argv[optind];

    struct qes_seq *seq = qes_seq_create();
    struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
    ssize_t res;
    while ((res = qes_seqfile_read(sf, seq)) > 0) {
        if (prefix != NULL || suffix != NULL) {
            seq_print_presuf(seq, stdout, prefix, prefixlen, suffix, suffixlen, fasta);
        }
    }

    qes_seqfile_destroy(sf);
    qes_seq_destroy(seq);
    return EXIT_SUCCESS;
}
