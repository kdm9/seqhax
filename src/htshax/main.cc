#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

struct prog {
    const char *name;
    int (*mainfunc)(int, char **);
    const char *desc;
};

int bcfhist_main(int argc, char *argv[]);


static const struct prog programs[] = {
    {"bcfhist",     bcfhist_main,       "Calculate histograms on MAF allele freq and missingness across samples"},
    {NULL,          NULL,               NULL}
};

void
htshax_usage()
{
    cerr << "Usage" << endl
         << "    htshax PROGRAM [options]" << endl
         << endl
         << "where PROGRAM is one of:" << endl
         << endl;;
    for (const struct prog *program = programs; program->name; program++) {
        cerr << program->name <<  " -- " << program->desc << endl;
    }
}


int
main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stdout, "htshax version: " SHX_VERSION "\n\n");
        htshax_usage();
        return EXIT_SUCCESS;
    }

    for (const struct prog *program = programs; program->name; program++) {
        if (strcmp(argv[1], program->name) == 0) {
            return program->mainfunc(argc - 1, argv + 1);
        }
    }

    /* If we get to here, no program matched the program name */
    fprintf(stderr, "Invalid program name '%s'\n", argv[1]);
    htshax_usage();
    return EXIT_FAILURE;
}

