/*
 * ============================================================================ *
 *       Filename:  unique_barcode.c
 *    Description:  
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
    struct qes_seq *seq = NULL;
    struct qes_seqfile *sf = NULL;
    ssize_t res = 0;
    size_t n_recs = 0;
    size_t len = 0;

    switch(getopt(argc, argv, "l:")) {
    case 'l':
        len = (size_t) strtoul(optarg, NULL, 10);
        break;
    case '?':
        fprintf(stderr, "Bad option %s\n\n", argv[optind]);
    default:
        fprintf(stderr, "Usage: %s -l <len> <fastqfile> \n", argv[0]);
        return EXIT_FAILURE;
    }
    if (argc <= optind) {
        fprintf(stderr, "Must provide input file\n\n");
        fprintf(stderr, "Usage: %s -l <len> <fastqfile> \n", argv[0]);
        return EXIT_FAILURE;
    }
    /* Make the prefix string */
    fprintf(stderr, "Removing '%zu' bases from each read of %s\n", len,
            argv[optind]);
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
        fprint_seq_without_prefix(stdout, seq, len);
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
