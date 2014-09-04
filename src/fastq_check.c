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
    struct qes_seq *seq = qes_seq_create();
    struct qes_seqfile *sf = NULL;
    ssize_t res = 0;
    size_t n_recs = 0;
    size_t seq_len = 0;
    int ok = 1;

    if (argc < 2) {
        fprintf(stderr, "Usage: fastq_check <fastqfile>\n");
        exit(127);
    }
    sf = qes_seqfile_create(argv[1], "r");
    while (res != EOF) {
        res = qes_seqfile_read(sf, seq);
        if (res < -2) {
            printf("[main] Error '%s' at read %zu\n", error_strings[-res],
                   ++n_recs);
            ok = 0;
            continue;
        } else if (res < 1) {
            break;
        }
        seq_len += res;
        if (++n_recs % 1000==0) {
            printf("[main] Processed %.3fM seqs\r", (float)n_recs/1000000.0);
        }
    }
    printf("[main] Processed %zu sequences\n", n_recs);
    printf("[main] Total sequence length is %zu base pairs\n", seq_len);
    if (ok) {
        printf("[main] No errors detected in %s\n", argv[1]);
    } else {
        printf("[main] Errors were detected in %s\n", argv[1]);
    }
    qes_seqfile_destroy(sf);
    qes_seq_destroy(seq);
    exit(ok ? 0 : 1);
}
