/*
 * ============================================================================
 *
 *       Filename:  fastq_check.c
 *
 *    Description:  Check a fastq for weird errors
 *
 *        Created:  04/09/14 13:51:28
 *        License:  GPLv3+
 *       Compiler:  gcc, clang
 *   Compile With:  cc -O3 -Wall -Wextra -std=c99 -lqes -o fastq_check fastq_check.c
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <qes_seqfile.h>

/* Index this with -res in the below loop */
static const char *error_strings[] = {
    "", /* 0 */
    "End of file!", /* 1, i.e. EOF */
    "Misc error in libqes", /* 2, i.e. misc error */
    "Bad header line",
    "Bad sequence line",
    "Bad qual header line",
    "Bad qual line",
    "Seq and qual are different lengths",
};

int
main(int argc, char *argv[])
{
    int all_ok = 1;
    int iii = 0;

    if (argc < 2) {
        fprintf(stderr, "Usage: fastq_check <fastqfile> [<fastqfile> [...]]\n");
        exit(127);
    }
    fprintf(stdout, "File\tSequences\tBasepairs\tOk?\n");
    for (iii = 1; iii < argc; iii++ ) {
        struct qes_seq *seq = qes_seq_create();
        struct qes_seqfile *sf = NULL;
        ssize_t res = 0;
        size_t n_recs = 0;
        size_t seq_len = 0;
        int ok = 1;
        int chars_to_wipe = 0;
        int jjj = 0;

        sf = qes_seqfile_create(argv[iii], "r");
        while (res != EOF) {
            res = qes_seqfile_read(sf, seq);
            if (res < -2) {
                fprintf(stderr, "Error '%s' at read %zu in %s\n",
                        error_strings[-res], ++n_recs, argv[iii]);
                ok = 0;
                all_ok = 0;
                continue;
            } else if (res < 1) {
                break;
            }
            seq_len += res;
            if (++n_recs % 100000==0) {
                chars_to_wipe = fprintf(stderr, "Processed %.1fM seqs of %s\r",
                                        (float)n_recs/1000000.0, argv[iii]);
                fflush(stderr);
            }
        }
        for (jjj = 0; jjj < chars_to_wipe; jjj++) {
            fprintf(stderr, " ");
        }
        fprintf(stderr, "\r");
        fflush(stderr);
        fprintf(stdout, "%s\t%zu\t%zu\t%s\n", argv[iii], n_recs, seq_len,
                ok ? "Yes" : "No");
        qes_seqfile_destroy(sf);
        qes_seq_destroy(seq);
    }
    exit(all_ok ? 0 : 1);
}
