/*
 * ============================================================================
 *
 *       Filename:  seqhax.c
 *
 *    Description:  seqhax
 *
 *        Created:  20/12/14 21:24:53
 *        License:  GPLv3+
 *       Compiler:  gcc, clang
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include "seqhax.h"

int
fprint_seq_with_prefix(FILE *stream, struct qes_seq *seq, char *prefix,
                       size_t len)
{
    size_t iii = 0;

    if (stream == NULL || seq == NULL || prefix == NULL) {
        return 1;
    }
    if (seq->qual.l) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.s, stream);
    if (seq->comment.s) {
        fputc(' ', stream);
        fputs(seq->comment.s, stream);
    }
    fputc('\n', stream);
    fputs(prefix, stream);
    fputs(seq->seq.s, stream);
    fputc('\n', stream);
    if (seq->qual.l) {
        fputs("+\n", stream);
        /* qual score of I is valid for all encodings */
        for (iii = 0; iii < len; iii ++) {
            fputc('I', stream);
        }
        fputs(seq->qual.s, stream);
        fputc('\n', stream);
    }
    return 0;
}

int
fprint_seq_without_prefix(FILE *stream, struct qes_seq *seq, size_t prefix_len)
{
    if (stream == NULL || seq == NULL) {
        return 1;
    }
    if (seq->qual.l) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.s, stream);
    if (seq->comment.s) {
        fputc(' ', stream);
        fputs(seq->comment.s, stream);
    }
    fputc('\n', stream);
    fputs(seq->seq.s + prefix_len, stream);
    fputc('\n', stream);
    if (seq->qual.l) {
        fputs("+\n", stream);
        fputs(seq->qual.s + prefix_len, stream);
        fputc('\n', stream);
    }
    return 0;
}

int
fprint_seq_with_suffix(FILE *stream, struct qes_seq *seq, char *suffix,
                       size_t len)
{
    size_t iii = 0;

    if (stream == NULL || seq == NULL || suffix == NULL) {
        return 1;
    }
    if (seq->qual.l) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.s, stream);
    if (seq->comment.s) {
        fputc(' ', stream);
        fputs(seq->comment.s, stream);
    }
    fputc('\n', stream);
    fputs(suffix, stream);
    fputs(seq->seq.s, stream);
    fputc('\n', stream);
    if (seq->qual.l) {
        fputs("+\n", stream);
        /* IIII is valid for all encodings */
        for (iii = 0; iii < len; iii ++) {
            fputc('I', stream);
        }
        fputs(seq->qual.s, stream);
    }
    fputc('\n', stream);
    return 0;
}

int
fprint_seq_without_suffix(FILE *stream, struct qes_seq *seq, size_t suffix_len)
{
    size_t trunc_idx = 0;
    char trunc_chr = '\0';

    if (stream == NULL || seq == NULL) {
        return 1;
    }
    trunc_idx = seq->seq.l - suffix_len - 1;
    if (seq->qual.l) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.s, stream);
    if (seq->comment.s) {
        fputc(' ', stream);
        fputs(seq->comment.s, stream);
    }
    fputc('\n', stream);
    trunc_chr = seq->seq.s[trunc_idx];
    seq->seq.s[trunc_idx] = '\0';
    fputs(seq->seq.s, stream);
    seq->seq.s[trunc_idx] = trunc_chr;
    fputc('\n', stream);
    if (seq->qual.l) {
        fputs("+\n", stream);
        trunc_chr = seq->qual.s[trunc_idx];
        seq->qual.s[trunc_idx] = '\0';
        fputs(seq->qual.s, stream);
        seq->qual.s[trunc_idx] = trunc_chr;
    }
    fputc('\n', stream);
    return 0;
}
