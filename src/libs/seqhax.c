#include "seqhax.h"

int
fprint_seq_with_prefix(FILE *stream, struct qes_seq *seq, char *prefix,
                       size_t len)
{
    size_t iii = 0;

    if (stream == NULL || seq == NULL || prefix == NULL) {
        return 1;
    }
    if (seq->qual.len) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.str, stream);
    if (seq->comment.str) {
        fputc(' ', stream);
        fputs(seq->comment.str, stream);
    }
    fputc('\n', stream);
    fputs(prefix, stream);
    fputs(seq->seq.str, stream);
    fputc('\n', stream);
    if (seq->qual.len) {
        fputs("+\n", stream);
        /* qual score of I is valid for all encodings */
        for (iii = 0; iii < len; iii ++) {
            fputc('I', stream);
        }
        fputs(seq->qual.str, stream);
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
    if (seq->qual.len) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.str, stream);
    if (seq->comment.str) {
        fputc(' ', stream);
        fputs(seq->comment.str, stream);
    }
    fputc('\n', stream);
    fputs(seq->seq.str + prefix_len, stream);
    fputc('\n', stream);
    if (seq->qual.len) {
        fputs("+\n", stream);
        fputs(seq->qual.str + prefix_len, stream);
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
    if (seq->qual.len) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.str, stream);
    if (seq->comment.str) {
        fputc(' ', stream);
        fputs(seq->comment.str, stream);
    }
    fputc('\n', stream);
    fputs(suffix, stream);
    fputs(seq->seq.str, stream);
    fputc('\n', stream);
    if (seq->qual.len) {
        fputs("+\n", stream);
        /* IIII is valid for all encodings */
        for (iii = 0; iii < len; iii ++) {
            fputc('I', stream);
        }
        fputs(seq->qual.str, stream);
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
    trunc_idx = seq->seq.len - suffix_len - 1;
    if (seq->qual.len) {
        fputc('@', stream);
    } else {
        fputc('>', stream);
    }
    fputs(seq->name.str, stream);
    if (seq->comment.str) {
        fputc(' ', stream);
        fputs(seq->comment.str, stream);
    }
    fputc('\n', stream);
    trunc_chr = seq->seq.str[trunc_idx];
    seq->seq.str[trunc_idx] = '\0';
    fputs(seq->seq.str, stream);
    seq->seq.str[trunc_idx] = trunc_chr;
    fputc('\n', stream);
    if (seq->qual.len) {
        fputs("+\n", stream);
        trunc_chr = seq->qual.str[trunc_idx];
        seq->qual.str[trunc_idx] = '\0';
        fputs(seq->qual.str, stream);
        seq->qual.str[trunc_idx] = trunc_chr;
    }
    fputc('\n', stream);
    return 0;
}
