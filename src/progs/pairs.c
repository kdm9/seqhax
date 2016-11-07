/*******************************************************************************
*                  pairs -- (De)interleave paired end reads                   *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
pairs_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax pairs [options] FILE [FILE2]\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -f         Force output to existing files.\n");
    fprintf(stream, "    -l LENGTH  Minimum length of each read. [default 1]\n");
    fprintf(stream, "    -1 FILE    Pair first mate output\n");
    fprintf(stream, "    -2 FILE    Pair second mate output\n");
    fprintf(stream, "    -u FILE    Unpaired read output\n");
    fprintf(stream, "    -s FILE    Strict interleaved output\n");
    fprintf(stream, "    -b FILE    Broken paired output\n");
    fprintf(stream, "\n");
    fprintf(stream, "Output files are NOT compressed. To apply compression, please use\n");
    fprintf(stream, "Subprocess streams, e.g.:\n");
    fprintf(stream, "\n");
    fprintf(stream, "  seqhax pairs -1 >(gzip > read1.fq.gz) -2 >(gzip > read2.fq.gz) \\\n");
    fprintf(stream, "      -u >(gzip > unapired.fq.gz) reads.fq.gz\n");
    fprintf(stream, "\n");
    fprintf(stream, "One can of course use other compression algorithms, e.g zstd.\n");
    fprintf(stream, "\n");
    fprintf(stream, "FILE should be a sequence file in FASTA or FASTQ format.\n");
    fprintf(stream, "To accept reads from standard input, use '/dev/stdin' as\n");
    fprintf(stream, "the input file. To output to standard output, use '-'.\n");
}
static const char *pairs_optstr = "fl:1:2:u:s:b:h";

/*******************************************************************************
*                                   Helpers                                   *
*******************************************************************************/


bool
ispaired(struct qes_seq *r1, struct qes_seq *r2)
{
    if (!qes_seq_ok(r1) || !qes_seq_ok(r2)) return false;
    puts("1");
    if (r1->name.len == 0 || r2->name.len == 0) return false;
    puts("2");
    if (r1->name.len != r2->name.len) return false;
    puts("3");
    for (size_t i = 0; i < r1->name.len; i++) {
        if (r1->name.str[i] != r2->name.str[i]) {
            puts("4");
            if (i == r1->name.len - 2 && r1->name.str[i] == '/') {
                puts("t");
                return true;
            }
            puts("5");
            return false;
        }
        puts("6");
        return true;
    }
    puts("7");
    return true;
}



/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
pairs_main(int argc, char *argv[])
{
    int c = 0;
    ssize_t min_len = 1;
    const char *infile1 = NULL;
    const char *infile2 = NULL;
    const char *ilfile = NULL;
    const char *r1file = NULL;
    const char *r2file = NULL;
    const char *rsfile = NULL;
    bool strictpaired = false;
    char outmode[2] = "x";

    /* Parse CLI options */
    while ((c = getopt(argc, argv, pairs_optstr)) >= 0) {
        switch (c) {
        case 'f':
            outmode[0] = 'w';
            break;
        case 'b':
            strictpaired = false;
            ilfile = optarg;
            break;
        case 's':
            strictpaired = true;
            ilfile = optarg;
            break;
        case '1':
            r1file = optarg;
            break;
        case '2':
            r2file = optarg;
            break;
        case 'u':
            rsfile = optarg;
            break;
        case 'l':
            min_len = (ssize_t) atol(optarg);
            break;
        case 'h':
            pairs_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            pairs_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (optind >= argc) {
        pairs_usage(stderr);
        return EXIT_FAILURE;
    }
    infile1 = argv[optind];
    if (optind + 1 < argc) {
        infile2 = argv[optind + 1];
    }

    /***************************************************************************
    *                               Open output files                          *
    ***************************************************************************/
    FILE *r1fp = NULL;
    FILE *r2fp = NULL;
    FILE *rsfp = NULL;
    if (ilfile != NULL) {
        FILE *ilfp;
        if (strcmp(ilfile, "-") == 0 || strcmp(ilfile, "/dev/stdout") == 0) {
            ilfp = stdout;
        } else {
            ilfp = fopen(ilfile, outmode);
        }
        if (ilfp == NULL) {
            fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f).\n",
                    ilfile);
            fputc('\n', stderr);
            pairs_usage(stderr);
            return EXIT_FAILURE;
        }
        r1fp = ilfp;
        r2fp = ilfp;
        rsfp = ilfp;
    } else {
        if (r1file == NULL || r2file == NULL || rsfile == NULL) {
            fputs("To split reads to separate files, -1, -2 and -s must all be specified.\n", stderr);
            fputs("Use /dev/null as filename to discard.\n", stderr);
            fputc('\n', stderr);
            pairs_usage(stderr);
            return EXIT_FAILURE;
        }

        r1fp = fopen(r1file, outmode);
        if (r1fp == NULL) {
            fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                    r1file);
        }
        r2fp = fopen(r2file, outmode);
        if (r2fp == NULL) {
            fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                    r2file);
        }
        rsfp = fopen(rsfile, outmode);
        if (rsfp == NULL) {
            fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                    rsfile);
        }
        if (r1fp == NULL || r2fp == NULL || rsfp == NULL) {
            if (r1fp != NULL) fclose(r1fp);
            if (r2fp != NULL) fclose(r2fp);
            if (rsfp != NULL) fclose(rsfp);
            fputc('\n', stderr);
            pairs_usage(stderr);
            return EXIT_FAILURE;
        }
    }

    /***************************************************************************
    *                         Open inputs, main loop                          *
    ***************************************************************************/

    struct qes_seq *r1 = qes_seq_create();
    struct qes_seq *r2 = qes_seq_create();
    struct qes_seqfile *sf1 = qes_seqfile_create(infile1, "r");
    struct qes_seqfile *sf2 = sf1;
    if (infile2 != NULL) {
        sf2 = qes_seqfile_create(infile2, "r");
    }
    int fasta = -1;
    ssize_t npair = 0, nr1s = 0, nr2s = 0, nsingle=0, nfail = 0;
    while (1) {
        ssize_t res1 = qes_seqfile_read(sf1, r1);
        ssize_t res2 = qes_seqfile_read(sf2, r2);
        if (res1 < 1 &&  res2 < 1) {
            break;
        }

        bool r1pass = res1 >= min_len;
        bool r2pass = res2 >= min_len;

        if (fasta < 0) fasta = r1->qual.len == 0;

        if (!r1pass && !r2pass) {
            nfail += 2;
            continue;
        }

        bool is_pair = ispaired(r1, r2);

        if (is_pair) {
            if (r1pass && r2pass) {
                qes_seq_print(r1, r1fp, fasta, 1);
                qes_seq_print(r2, r2fp, fasta, 2);
                npair ++;
            } else {
                if (strictpaired) {
                    if (!r1pass) {
                        qes_seq_fill_seq(r1, "N", 1);
                        qes_seq_fill_qual(r1, "I", 1);
                        nr2s++;
                    } else {
                        // r2 must have failed
                        qes_seq_fill_seq(r2, "N", 1);
                        qes_seq_fill_qual(r2, "I", 1);
                        nr1s++;
                    }
                    qes_seq_print(r1, r1fp, fasta, 1);
                    qes_seq_print(r2, r2fp, fasta, 2);
                } else {
                    if (r1pass) {
                        qes_seq_print(r1, rsfp, fasta, 1);
                        nr1s++;
                    }
                    if (r2pass) {
                        qes_seq_print(r2, rsfp, fasta, 2);
                        nr2s++;
                    }
                }
            }
        } else { // Not a pair
            if (r1pass) {
                qes_seq_print(r1, rsfp, fasta, 0);
                nsingle++;
                if (strictpaired) {
                    qes_seq_fill_seq(r1, "N", 1);
                    qes_seq_fill_qual(r1, "I", 1);
                    qes_seq_print(r1, rsfp, fasta, 0);
                }
            } else {
                nfail++;
            }
            if (r2pass) {
                qes_seq_print(r2, rsfp, fasta, 0);
                nsingle++;
                if (strictpaired) {
                    qes_seq_fill_seq(r2, "N", 1);
                    qes_seq_fill_qual(r2, "I", 1);
                    qes_seq_print(r2, rsfp, fasta, 0);
                }
            } else {
                nfail++;
            }
        }
    }
    fprintf(stderr, "---\nstats:\n");
    fprintf(stderr, "  proper_pairs: %zu\n", npair);
    fprintf(stderr, "  single_r1: %zu\n", nr1s);
    fprintf(stderr, "  single_r2: %zu\n", nr2s);
    fprintf(stderr, "  unpaired_reads: %zu\n", nsingle);
    fprintf(stderr, "  failed: %zu\n", nfail);

    if (sf1 != sf2) {
        qes_seqfile_destroy(sf2);
    }
    qes_seqfile_destroy(sf1);
    qes_seq_destroy(r1);
    qes_seq_destroy(r2);
    return EXIT_SUCCESS;
} 
