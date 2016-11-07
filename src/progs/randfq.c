/*******************************************************************************
*                  randfq -- Creates a random sequence file                   *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>

#include "qes_seq.h"
#include "pcg_variants.h"

static void
randfq_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax randfq [options]\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -n READS   Number of reads. Use 0 for infinite. [default 1000]\n");
    fprintf(stream, "    -s SEED    Seed for RNG.\n");
    fprintf(stream, "    -l LENGTH  Length of each read. [default 100]\n");
    fprintf(stream, "    -p         Paired reads [default false]\n");
    fprintf(stream, "    -f         Output as fasta (no qualities)\n");
}
static const char *randfq_optstr = "n:s:l:pfh";

static inline void
randseq(pcg32_random_t *rng, char *seq, size_t len)
{
    assert(rng != NULL);
    assert(seq != NULL);
    const char NT[] = "ACGT";
    uint32_t rand;

    for (size_t i = 0; i < len;) {
        rand = pcg32_random_r(rng);
        for (size_t shift = 0; shift < 16 && i < len; shift++) {
            seq[i++] = NT[rand & 0x3];
            rand >>= 2;
        }
    }
}

int
randfq_main(int argc, char *argv[])
{
    struct qes_seq *seq;
    size_t n_recs = 0;
    uint32_t seqlen = 100;
    uint64_t num_seqs = 1000;
    pcg32_random_t rng;
    uint64_t seed = 0;
    bool paired = false;
    bool fasta = false;
    int c;

    while ((c = getopt(argc, argv, randfq_optstr)) > 0) {
        switch (c) {
            case 'n':
                num_seqs = strtoull(optarg, NULL, 10);
                break;
            case 's':
                seed = strtoul(optarg, NULL, 10);
                if (errno != 0) {
                    fprintf(stderr, "Invalid seed '%s'\n\n", optarg);
                    randfq_usage(stderr);
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                seqlen = strtoull(optarg, NULL, 10);
                if (seqlen > UINT32_MAX) {
                    fprintf(stderr, "Invalid read length '%s'\n\n", optarg);
                    randfq_usage(stderr);
                    return EXIT_FAILURE;
                }
                break;
            case 'p':
                paired = true;
                break;
            case 'f':
                fasta = true;
                break;
            case 'h':
            default:
                randfq_usage(stderr);
                return EXIT_FAILURE;
                break;
        }
    }
    if (seed == 0) {
        seed = (uint64_t) time(NULL);
        // x-or with address in case multiple processes start at once;
        seed ^= (uint64_t) &seed;
    }
    pcg32_srandom_r(&rng, seed, seed + 12);

    // Init sequence
    char *sequence = calloc(seqlen + 1, 1);
    assert(sequence != NULL);

    seq = qes_seq_create();
    if (fasta) {
        qes_seq_fill_qual(seq, "", 0);
    } else {
        char *qual = calloc(seqlen + 1, 1);
        assert(qual != NULL);
        for (size_t i = 0; i < seqlen; i++) {
            qual[i] = 'I';
        }
        qes_seq_fill_qual(seq, qual, seqlen);
        free(qual);
    }

    for (; num_seqs == 0 || n_recs < num_seqs; n_recs++) {
        char readid[101] = "";
        int tag = 0;
        const size_t namelen = snprintf(readid, 100, "read:%zu", n_recs);
        qes_seq_fill_name(seq, readid, namelen);

        randseq(&rng, sequence, seqlen);
        qes_seq_fill_seq(seq, sequence, seqlen);
        
        if (paired) {
            tag = (n_recs % 2) + 1;
        }
        qes_seq_print(seq, stdout, fasta, tag);
    }
    qes_seq_destroy(seq);
    free(sequence);
    return EXIT_SUCCESS;
}
