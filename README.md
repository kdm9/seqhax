Seqhax
======

A collection of (mostly C) programs that are vaguely useful for sequence
analysis. None deserve their own git repo (yet). They may be useful to others,
or maybe not.

Installation
============

Install zlib version >=1.2.5, then:

    git clone https://github.com/kdmurray91/seqhax.git
    cd seqhax
    mkdir build && cd build
    cmake ..
    make
    make install

Any other issues, file a bug report on github, as it probably is one.

Documentation
=============

fastq-check
-----------

Check a fastq file for syntax errors. Also counts the number of records and
total sequence length in base-pairs.

Usage:

    fastq-check <fastq-file>  [<fastqfile> [...]]

It is a tab-seperated table with files as rows. I use it like this:

    fastq-check reads/*.fastq.gz | tee read_summary.tab | column -t

This saves a tsv file for further use/plotting in R or whatever, and prints a
nicely formatted table to the terminal. There's also a progress line printed to
stderr, so you know how far through the files it is.

Note that the syntax is checks is fairly basic: there must be four lines, with
header, sequence, a line that matches '\+.*\n', and the quality scores. That's
it. Technically valid fastqs that have multiline seq/quality or that have
multibyte quality scores etc are counted as invalid, as libqes doesn't handle
these for speed. Likewise, libqes doesn't do proper checking of valid quality
scores, basically any ASCII is counted as valid.

License
=======

GPL v3 (see ./LICENSE).

Copyright (c) 2014 Kevin Murray.

- `src/ssw.[ch]` is Copyright (c) 2012-1015 Boston College, and are part of
  [ssw](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library).
