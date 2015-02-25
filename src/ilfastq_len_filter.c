/*
 * ============================================================================
 *
 *       Filename:  ilfastq_len_filter.c
 *    Description:  Check a fastq for weird errors
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <qes_seqfile.h>
#include <qes_seq.h>

int
main(int argc, char *argv[])
{
    struct qes_seq *r1 = NULL;
    struct qes_seq *r2 = NULL;
    struct qes_seqfile *sf = NULL;
    size_t min_len = 0;
    size_t n_recs = 0;
    ssize_t res1 = 0;
    ssize_t res2 = 0;

    /* Parse CLI options */
    switch(getopt(argc, argv, "l:")) {
    case 'l':
        min_len = (size_t) atol(optarg);
        break;
    case '?':
        fprintf(stderr, "Bad option %s\n\n", argv[optind]);
    default:
        fprintf(stderr, "USAGE: %s -l <min_len> <ilfq_file>\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (argc <= optind) {
        fprintf(stderr, "Error: Must provide input file\n\n");
        fprintf(stderr, "USAGE: %s -l <min_len> <ilfq_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    r1 = qes_seq_create();
    r2 = qes_seq_create();
    sf = qes_seqfile_create(argv[optind], "r");
    while ((res1 = qes_seqfile_read(sf, r1)) > 0 &&
           (res2 = qes_seqfile_read(sf, r2)) > 0) {
        if (r1->seq.len >= min_len && r2->seq.len >= min_len) {
            qes_seq_print(r1, stdout);
            qes_seq_print(r2, stdout);
        }
        if (++n_recs % 100000==0) {
            fprintf(stderr, "Processed %.1fM seqs of %s\r",
                    (float)n_recs/1000000.0, argv[optind]);
            fflush(stderr);
        }
    }

    qes_seqfile_destroy(sf);
    qes_seq_destroy(r1);
    qes_seq_destroy(r2);
    return EXIT_FAILURE;
}
