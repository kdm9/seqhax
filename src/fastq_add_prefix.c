/*
 * ============================================================================ *
 *       Filename:  fastq_add_prefix.c
 *    Description:  Add a DNA prefix to each read
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */
#include <stdlib.h>

#include <qes.h>
#include <seqhax.h>
#include <getopt.h>

int
main(int argc, char *argv[])
{
    char *prefix = NULL;
    size_t prefix_len = 0;
    struct qes_seq *seq = NULL;
    struct qes_seqfile *sf = NULL;
    ssize_t res = 0;
    size_t n_recs = 0;

    switch(getopt(argc, argv, "b:")) {
    case 'p':
        prefix = strdup(optarg);
        prefix_len = strlen(prefix);
        break;
    case '?':
        fprintf(stderr, "Bad option %s\n\n", argv[optind]);
        fprintf(stderr, "Usage: %s -p <prefix> <fastqfile> \n", argv[0]);
    default:
        fprintf(stderr, "Usage: %s -p <prefix> <fastqfile> \n", argv[0]);
        return 1;
    }
    if (argc <= optind) {
        fprintf(stderr, "Error: Must provide input file\n\n");
        fprintf(stderr, "Usage: %s -p <prefix> <fastqfile> \n", argv[0]);
        return 1;
    }
    fprintf(stderr, "Adding prefix '%s' to %s\n", prefix, argv[optind]);
    seq = qes_seq_create();
    sf = qes_seqfile_create(argv[optind], "r");
    while (res != EOF) {
        res = qes_seqfile_read(sf, seq);
        if (res == EOF) {
            break;
        } else if (res < 1) {
            fprintf(stderr, "Broken file at seq %zu\n", n_recs);
            break;
        }
        fprint_seq_with_prefix(stdout, seq, prefix, prefix_len);
        if (++n_recs % 100000==0) {
            fprintf(stderr, "Processed %.1fM seqs of %s\r",
                                    (float)n_recs/1000000.0, argv[optind]);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\nDone!\n");
    qes_seq_destroy(seq);
    qes_seqfile_destroy(sf);
    return 0;
}
