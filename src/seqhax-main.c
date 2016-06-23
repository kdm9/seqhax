#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "progs/progs.h"

static const char *programs[] = {
    "randfq  -- generates random sequence files",
    NULL
};

struct seqhax_prog {
    const char *name;
    int (*mainfunc)(int, char **);
};

static const struct seqhax_prog program_mains[] = {
    {"randfq", randfq_main},
    {NULL, NULL}
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
