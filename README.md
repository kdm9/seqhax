Seqhax
======

A seqtk-style toolkit for sequence analysis. By no means feature complete. In
fact largely contains features other authors have not implemented in their
respective tools.

Installation
============

Install zlib version >=1.2.5, then:

```bash
git clone https://github.com/kdmurray91/seqhax.git
cd seqhax
mkdir build && cd build
cmake ..
make
make install
```

To make static binaries, one can use `cmake -DSTATIC_BUILD=On ..` in the above series of commands.

Any other issues, file a [bug report on github](https://github.com/kdmurray91/seqhax/issues).

Documentation
=============

There are two commands in this package: `seqhax` and `htshax`. 

## `seqhax`

The `seqhax` command contains various sequence file manipulation subcommands. The commands, along with a synopsis
of their actions, are displayed when one types `seqhax` with no arguments. 
At the time of writing, these were:

```
$ seqhax
USAGE:
    seqhax PROGRAM [options]

where PROGRAM is one of:
    anon       -- Rename sequences by a sequential number
    convert    -- Convert between FASTA and FASTQ formats
    clihist    -- Records a histogram of stdin.
    filter     -- Filter reads from a sequence file
    pairs      -- (De)interleave paired end reads
    pecheck    -- Check that paired end reads match properly
    preapp     -- Prepend or append string to sequences
    randseq    -- Generate a random sequence file
    rebarcode  -- Add index sequences from header to the start of reads
    stats      -- Basic statistics about sequence files
    trunc      -- Truncate sequences
```

The usage of each subcommand can be obtained using the `-h` flag to that
subcommand, e.g. `seqhax preapp -h`. 

Some brief description of the more involved commands is below.

#### `pairs`

Interleaves or de-interleaves paired reads. Converts between the following
forms:

- Separate R1/R2 paired files and single read read file.
- "Strict" interleaved file, where failed/missing reads are replaced with a
  single 'N'.
- "Broken paired" interleaved files, where failed/missing reads are simply not
  present.

#### `pecheck`

Checks that read pairs are correctly matched, between split (R1 & R2) files, or
interleaved files. Optionally, can be used to join multiple R1/R2 from the
sample sample into a single interleaved file, while checking read IDs match.

## `htshax`

`htshax` consists of various utilities for the handling of HTS formats via HTSlib. These are mostly things that I wish samtools/bcftools would do, but don't.

```
$ htshax
Usage
    htshax PROGRAM [options]

where PROGRAM is one of:

bcfhist -- Calculate histograms on MAF allele freq and missingness across samples
```

#### `bcfhist`

Create histograms over various V/BCF metrics used in variant filtering (`QUAL`, `ALLELE_FREQ`, `F_MISSING`, `INFO/DP`). Nominally bcftools stats could do this, but it doesn't summarise missingness, so I wrote this quick tool to do so.

# License

GPL v3 (see ./LICENSE).

Copyright (c) 2014-2023 Kevin Murray.
