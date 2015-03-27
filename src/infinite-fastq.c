/*
 * ============================================================================
 *
 *       Filename:  infinite-fastq.c
 *    Description:  Make an infinite FASTQ file
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <qes_seq.h>

int
main(int argc, char *argv[])
{
    struct qes_seq *seq;
    size_t n_recs = 0;
    char num[22] = "";
    size_t len = 0;

    seq = qes_seq_create();
    qes_seq_fill(seq, "", "", "ACGTACGTACGTACGTACGT", "22222222222222222222");
    while (1){
        len = snprintf(num, 22, "%zu", n_recs++);
        qes_seq_fill_name(seq, num, len);
        qes_seq_fill_comment(seq, num, len);
        qes_seq_print(seq, stdout);
    }
    qes_seq_destroy(seq);
    return EXIT_SUCCESS;
}
