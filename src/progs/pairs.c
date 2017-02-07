/*******************************************************************************
*                  pairs -- (De)interleave paired end reads                   *
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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
    fprintf(stream, "    -p FILE    Interleaved pairs-only output\n");
    fprintf(stream, "    -u FILE    Unpaired read output\n");
    fprintf(stream, "    -s FILE    \"Strict interleaved\" output, all reads\n");
    fprintf(stream, "    -b FILE    \"Broken paired\" output, all reads\n");
    fprintf(stream, "    -y FILE    Output statistics to FILE.\n");
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
static const char *pairs_optstr = "fl:1:2:p:u:s:b:y:h";

/*******************************************************************************
*                                   Helpers                                   *
*******************************************************************************/

bool
isfilepipe(const char *filename)
{
    struct stat res;
    if (stat(filename, &res) != 0) return 0;
    return (S_ISFIFO(res.st_mode));
}

bool
ispaired(struct qes_seq *r1, struct qes_seq *r2)
{
    if (!qes_seq_ok(r1) || !qes_seq_ok(r2)) return false;
    if (r1->name.len == 0 || r2->name.len == 0) return false;
    for (size_t i = 0; i < r1->name.len; i++) {
        if (r1->name.str[i] != r2->name.str[i]) {
            if (i > 1 && r1->name.str[i-1] == '/' && r2->name.str[i-1] == '/'
                      && r1->name.str[i] == '1' && r2->name.str[i]=='2') {
                // Old illumina style:
                // @HWIblah/1 and @HWIblah/2
                return true;
            }
            return false;
        }
    }
    // New illumina style:
    // @HWIblah 1:blah and @HWIblah 2:blah
    // Name str stops at first space char. Names identical == pairs
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
    const char *paironlyfile = NULL;
    const char *statsfile = NULL;
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
        case 'p':
            paironlyfile = optarg;
            break;
        case 'y':
            statsfile = optarg;
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
        } else  {
            const char *om = isfilepipe(ilfile) ? "w" : outmode;
            ilfp = fopen(ilfile, om);
        }
        if (ilfp == NULL) {
            fprintf(stderr, "Could not open '%s' for IL output. Perhaps it already exists? (use -f).\n",
                    ilfile);
            fputc('\n', stderr);
            pairs_usage(stderr);
            return EXIT_FAILURE;
        }
        r1fp = ilfp;
        r2fp = ilfp;
        rsfp = ilfp;
    } else {
        if ((paironlyfile == NULL && (r1file == NULL || r2file == NULL))) {
            fputs("Either -s, -b, -p, or  -1 AND -2 must be specified.\n", stderr);
            fputs("Use /dev/null as filename to discard.\n", stderr);
            fputc('\n', stderr);
            pairs_usage(stderr);
            return EXIT_FAILURE;
        }

        if (paironlyfile != NULL) {
            const char *om = isfilepipe(paironlyfile) ? "w" : outmode;
            r1fp = fopen(paironlyfile, om);
            if (r1fp == NULL) {
                fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                        paironlyfile);
            }
            r2fp = r1fp;
        } else {
            const char *om = isfilepipe(r1file) ? "w" : outmode;
            r1fp = fopen(r1file, om);
            if (r1fp == NULL) {
                fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                        r1file);
            }
            om = isfilepipe(r2file) ? "w" : outmode;
            r2fp = fopen(r2file, om);
            if (r2fp == NULL) {
                fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                        r2file);
            }
        }
        if (rsfile != NULL) {
            const char *om = isfilepipe(rsfile) ? "w" : outmode;
            rsfp = fopen(rsfile, om);
            if (rsfp == NULL) {
                fprintf(stderr, "Could not open '%s' for output. Perhaps it already exists? (use -f)\n",
                        rsfile);
            }
        }
        if (r1fp == NULL || r2fp == NULL) {
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
    bool seekpair = false;
    size_t nread = 0;
    while (1) {
        ssize_t res1 = 0;
        if (!seekpair) {
            res1 = qes_seqfile_read(sf1, r1);
        }
        ssize_t res2 = qes_seqfile_read(sf2, r2);
        if (seekpair && res2 < 1) {
            qes_seq_print(r1, rsfp, fasta, 0);
            nsingle++;
            if (strictpaired) {
                qes_seq_fill_seq(r1, "N", 1);
                qes_seq_fill_qual(r1, "I", 1);
                qes_seq_print(r1, rsfp, fasta, 0);
            }
            break;
        }
        if (res1 < 1 &&  res2 < 1) {
            break;
        }
        bool r1pass = seekpair || res1 >= min_len;
        bool r2pass = res2 >= min_len;
        nread += seekpair ? 1 : 2;
        seekpair = false;

        // Detect if we should use fasta from first 2 read qualities
        if (fasta < 0) fasta = r1->qual.len == 0 && r2->qual.len == 0;

        bool is_pair = ispaired(r1, r2);
        if (!is_pair) {
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
                // Swap so R2 is now R1
                struct qes_seq *tmpread = r1;
                r1 = r2;
                r2 = tmpread;
                seekpair = true;
            } else {
                nfail++;
            }
            continue;
        }


        if (!r1pass && !r2pass) {
            nfail += 2;
            continue;
        }


        if (r1pass && r2pass) {
            qes_seq_print(r1, r1fp, fasta, 1);
            qes_seq_print(r2, r2fp, fasta, 2);
            npair++;
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
    }
    fprintf(stderr, "Finished, processed %zu reads\n", nread);
    FILE *statsfp = stderr;
    if (statsfile != NULL) {
        statsfp = fopen(statsfile, "w");
        if (statsfp == NULL) {
            fprintf(stderr, "ERROR: Could not open stats file '%s'. "
                    "Using stderr.\n", statsfile);
            statsfp = stderr;
        }
    }
    fprintf(statsfp, "---\n");
    fprintf(statsfp, "files:\n");
    fprintf(statsfp, "  - \"%s\"\n", infile1);
    if (infile2) {
        fprintf(statsfp, "  - \"%s\"\n", infile2);
    }
    fprintf(statsfp, "stats:\n");
    fprintf(statsfp, "  total_reads: %zu\n", nread);
    fprintf(statsfp, "  proper_pairs: %zu\n", npair);
    fprintf(statsfp, "  single_r1: %zu\n", nr1s);
    fprintf(statsfp, "  single_r2: %zu\n", nr2s);
    fprintf(statsfp, "  unpaired_reads: %zu\n", nsingle);
    fprintf(statsfp, "  failed: %zu\n", nfail);

    if (sf1 != sf2) {
        qes_seqfile_destroy(sf2);
    }
    qes_seqfile_destroy(sf1);
    qes_seq_destroy(r1);
    qes_seq_destroy(r2);
    fcloseall();
    return EXIT_SUCCESS;
}
