#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *programs[] = {
    "anon       -- Rename sequences by a sequential number",
    "convert    -- Convert between FASTA and FASTQ formats",
    /* seqhax cat is unadvertised, since it's mainly for testing */
    "clihist    -- Records a histogram of stdin.",
    "filter     -- Filter reads from a sequence file",
    "pairs      -- (De)interleave paired end reads",
    "pecheck    -- Check that paired end reads match properly",
    "preapp     -- Prepend or append string to sequences",
    "randseq    -- Generate a random sequence file",
    "rebarcode  -- Add index sequences from header to the start of reads",
    "stats      -- Basic statistics about sequence files",
    "trunc      -- Truncate sequences",
    NULL
};

struct seqhax_prog {
    const char *name;
    int (*mainfunc)(int, char **);
};

extern "C" {
int anon_main(int argc, char *argv[]);
int convert_main(int argc, char *argv[]);
int filter_main(int argc, char *argv[]);
int pairs_main(int argc, char *argv[]);
int preapp_main(int argc, char *argv[]);
int randseq_main(int argc, char *argv[]);
int trunc_main(int argc, char *argv[]);
}

int cat_main(int argc, char *argv[]);
int clihist_main(int argc, char *argv[]);
int pecheck_main(int argc, char *argv[]);
int rebarcode_main(int argc, char *argv[]);
int stats_main(int argc, char *argv[]);
int grepbed_main(int argc, char *argv[]);


static const struct seqhax_prog program_mains[] = {
    {"anon",            anon_main},
    {"cat",             cat_main},
    {"convert",         convert_main},
    {"filter",          filter_main},
    {"pairs",           pairs_main},
    {"preapp",          preapp_main},
    {"randseq",         randseq_main},
    {"stats",           stats_main},
    {"trunc",           trunc_main},
    {"rebarcode",       rebarcode_main},
    {"pecheck",         pecheck_main},
    {"clihist",         clihist_main},
    {"grepbed",         grepbed_main},
    {NULL,              NULL}
};


void
seqhax_usage(FILE *stream)
{
    const char **progname = programs;
    fprintf(stream, "USAGE:\n");
    fprintf(stream, "    seqhax PROGRAM [options]\n");
    fprintf(stream, "\n");
    fprintf(stream, "where PROGRAM is one of:\n");
    while (*progname) {
        fprintf(stream, "    %s\n", *progname);
        progname++;
    }
}


int
main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stdout, "seqhax version: " SHX_VERSION "\n\n");
        seqhax_usage(stdout);
        return EXIT_SUCCESS;
    }

    const struct seqhax_prog *program = program_mains;
    for (program = program_mains; program->name; program++) {
        if (strcmp(argv[1], program->name) == 0) {
            return program->mainfunc(argc - 1, argv + 1);
        }
    }

    /* If we get to here, no program matched the program name */
    fprintf(stderr, "Invalid program name '%s'\n", argv[1]);
    seqhax_usage(stderr);
    return EXIT_FAILURE;
}
