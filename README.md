Seqhax
======

A collection of (mostly C) programs that are vaguely useful for sequence
analysis. None deserve their own git repo (yet). They may be useful to others,
or maybe not.

Installation
============

Install zlib and [libqes](https://github.com/kdmurray91/libqes).

    git clone https://github.com/kdmurray91/seqhax.git
    cd seqhax
    mkdir build && cd build
    cmake ..
    make
    make install

Documentation
=============

fastq_check
-----------

Check a fastq file for syntax errors. Also counts the number of records and
total sequence length in base-pairs.

Usage:

    fastq_check <fastq_file>

The output is relatively self explanatory :)


License
=======

GPL v3 (see ./LICENSE).

Copyright (c) 2014 Kevin Murray.

- `src/ssw.[ch]` is Copyright (c) 2012-1015 Boston College, and are part of
  [ssw](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library).
- `src/ksw.[ch]`, `src/kseq.h` and `src/trimadap.c` are Copyright (c) Heng Li
  2008-2014, with modifications by Kevin Murray, and were part of
  [seqtk](https://github.com/lh3/seqtk).
