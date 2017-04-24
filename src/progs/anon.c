/*******************************************************************************
*             anon -- Anonymise sequences by renaming 1, 2, etc.              *
*******************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>

#include "qes_seqfile.h"


static void
anon_usage(FILE *stream)
{
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax anon [options] FILE\n");
    fprintf(stream, "\n");
    fprintf(stream, "OPTIONS:\n");
    fprintf(stream, "    -x     Use base-16 sequence IDs.\n");
    fprintf(stream, "    -p     Treat reads as pairs, add /1 or /2 to headers.\n");
}
static const char *anon_optstr = "xp";


#define NAMEBUFLEN 100

/*******************************************************************************
*                                    MAIN                                     *
*******************************************************************************/
int
anon_main(int argc, char *argv[])
{
    int c = 0;
    const char *filename = NULL;
    bool hex = false, paired = false;

    /* Parse CLI options */
    while ((c = getopt(argc, argv, anon_optstr)) >= 0) {
        switch (c) {
        case 'p':
            paired = true;
            break;
        case 'x':
            hex = true;
            break;
        case 'h':
            anon_usage(stdout);
            return EXIT_SUCCESS;
            break;
        case '?':
            anon_usage(stderr);
            return EXIT_FAILURE;
            break;
        }
    }
    if (optind >= argc) {
        anon_usage(stderr); return EXIT_FAILURE;
    }
    filename = argv[optind];

    struct qes_seq *seq = qes_seq_create();
    struct qes_seqfile *sf = qes_seqfile_create(filename, "r");
    ssize_t res;
    size_t readno = 0;
    int readinpair = 0;
    char namebuf[NAMEBUFLEN] = "";
    char commentbuf[] = "";
    while ((res = qes_seqfile_read(sf, seq)) > 0) {
        if (paired) {
            if (readinpair == 0) {
                readno++;
            }
        } else {
            readno++;
        }
        if (hex) {
            snprintf(namebuf, NAMEBUFLEN, "%zx", readno);
        } else {
            snprintf(namebuf, NAMEBUFLEN, "%zu", readno);
        }
        qes_seq_fill_name(seq, namebuf, NAMEBUFLEN);
        qes_seq_fill_comment(seq, commentbuf, 1);
        qes_seq_print(seq, stdout, !qes_seq_has_qual(seq), paired ? readinpair + 1 : 0);
        readinpair = (readinpair + 1) % 2;
    }

    qes_seqfile_destroy(sf);
    qes_seq_destroy(seq);
    return EXIT_SUCCESS;
}
