/*
 * ============================================================================ *
 *       Filename:  unique_barcode.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  17/09/14 14:49:11
 *       Revision:  none
 *        License:  GPLv3+
 *       Compiler:  gcc, clang
 *   Compile With:  cc -O3 -Wall -Wextra -std=c99 -o  unique_barcode.c
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */
#include <stdlib.h>

#include <qes.h>
#include <getopt.h>

static const char *reverse_enc[256] = {
    "AAAA", "CAAA", "TAAA", "GAAA", "ACAA", "CCAA", "TCAA", "GCAA",
    "ATAA", "CTAA", "TTAA", "GTAA", "AGAA", "CGAA", "TGAA", "GGAA",
    "AACA", "CACA", "TACA", "GACA", "ACCA", "CCCA", "TCCA", "GCCA",
    "ATCA", "CTCA", "TTCA", "GTCA", "AGCA", "CGCA", "TGCA", "GGCA",
    "AATA", "CATA", "TATA", "GATA", "ACTA", "CCTA", "TCTA", "GCTA",
    "ATTA", "CTTA", "TTTA", "GTTA", "AGTA", "CGTA", "TGTA", "GGTA",
    "AAGA", "CAGA", "TAGA", "GAGA", "ACGA", "CCGA", "TCGA", "GCGA",
    "ATGA", "CTGA", "TTGA", "GTGA", "AGGA", "CGGA", "TGGA", "GGGA",
    "AAAC", "CAAC", "TAAC", "GAAC", "ACAC", "CCAC", "TCAC", "GCAC",
    "ATAC", "CTAC", "TTAC", "GTAC", "AGAC", "CGAC", "TGAC", "GGAC",
    "AACC", "CACC", "TACC", "GACC", "ACCC", "CCCC", "TCCC", "GCCC",
    "ATCC", "CTCC", "TTCC", "GTCC", "AGCC", "CGCC", "TGCC", "GGCC",
    "AATC", "CATC", "TATC", "GATC", "ACTC", "CCTC", "TCTC", "GCTC",
    "ATTC", "CTTC", "TTTC", "GTTC", "AGTC", "CGTC", "TGTC", "GGTC",
    "AAGC", "CAGC", "TAGC", "GAGC", "ACGC", "CCGC", "TCGC", "GCGC",
    "ATGC", "CTGC", "TTGC", "GTGC", "AGGC", "CGGC", "TGGC", "GGGC",
    "AAAT", "CAAT", "TAAT", "GAAT", "ACAT", "CCAT", "TCAT", "GCAT",
    "ATAT", "CTAT", "TTAT", "GTAT", "AGAT", "CGAT", "TGAT", "GGAT",
    "AACT", "CACT", "TACT", "GACT", "ACCT", "CCCT", "TCCT", "GCCT",
    "ATCT", "CTCT", "TTCT", "GTCT", "AGCT", "CGCT", "TGCT", "GGCT",
    "AATT", "CATT", "TATT", "GATT", "ACTT", "CCTT", "TCTT", "GCTT",
    "ATTT", "CTTT", "TTTT", "GTTT", "AGTT", "CGTT", "TGTT", "GGTT",
    "AAGT", "CAGT", "TAGT", "GAGT", "ACGT", "CCGT", "TCGT", "GCGT",
    "ATGT", "CTGT", "TTGT", "GTGT", "AGGT", "CGGT", "TGGT", "GGGT",
    "AAAG", "CAAG", "TAAG", "GAAG", "ACAG", "CCAG", "TCAG", "GCAG",
    "ATAG", "CTAG", "TTAG", "GTAG", "AGAG", "CGAG", "TGAG", "GGAG",
    "AACG", "CACG", "TACG", "GACG", "ACCG", "CCCG", "TCCG", "GCCG",
    "ATCG", "CTCG", "TTCG", "GTCG", "AGCG", "CGCG", "TGCG", "GGCG",
    "AATG", "CATG", "TATG", "GATG", "ACTG", "CCTG", "TCTG", "GCTG",
    "ATTG", "CTTG", "TTTG", "GTTG", "AGTG", "CGTG", "TGTG", "GGTG",
    "AAGG", "CAGG", "TAGG", "GAGG", "ACGG", "CCGG", "TCGG", "GCGG",
    "ATGG", "CTGG", "TTGG", "GTGG", "AGGG", "CGGG", "TGGG", "GGGG",
};

static int
fprint_seq_with_prefix(FILE *stream, struct qes_seq *seq, char *prefix)
{
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
        /* IIII is valid for all encodings */
        fputs("IIIIIIIIIIIIIIII", stream);
        fputs(seq->qual.s, stream);
    }
    fputc('\n', stream);
    return 0;
}

int
main(int argc, char *argv[])
{
    char prefix[17] = "";
    uint32_t file_id = 0;
    struct qes_seq *seq = NULL;
    struct qes_seqfile *sf = NULL;
    ssize_t res = 0;
    size_t n_recs = 0;
    int iii = 0;

    switch(getopt(argc, argv, "b:")) {
    case 'b':
        file_id = (uint32_t) strtoul(optarg, NULL, 10);
        break;
    case '?':
        fprintf(stderr, "Bad option %s\n\n", argv[optind]);
    default:
        fprintf(stderr, "Usage: unique_barcode -b <id> <fastqfile> \n");
        return 1;
    }
    if (argc <= optind) {
        fprintf(stderr, "Must provide input file\n\n");
        fprintf(stderr, "Usage: unique_barcode -b <id> <fastqfile> \n");
        return 1;
    }
    /* Make the prefix string */
    for (iii = 0; iii < 4; iii++) {
        memcpy(prefix + iii * 4, reverse_enc[(file_id & (0xff << 8 * iii)) >> iii * 8], 4);
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
        fprint_seq_with_prefix(stdout, seq, prefix);
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
