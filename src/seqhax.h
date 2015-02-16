/*
 * ============================================================================
 *
 *       Filename:  seqhax.h
 *    Description:  seqhax utilities
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>

#include <qes.h>

int fprint_seq_with_prefix     (FILE                   *stream,
                                struct qes_seq         *seq,
                                char                   *prefix,
                                size_t                  len);
int fprint_seq_with_suffix     (FILE                   *stream,
                                struct qes_seq         *seq,
                                char                   *suffix,
                                size_t                  len);
int fprint_seq_without_prefix  (FILE                   *stream,
                               struct qes_seq          *seq,
                               size_t                   prefix_len);
int fprint_seq_without_suffix (FILE                    *stream,
                               struct qes_seq          *seq,
                               size_t                   suffix_len);
