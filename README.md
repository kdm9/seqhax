Seqhax
======

A seqtk-style toolkit for sequence analysis. By no means feature complete. In
fact largely contains features other authors have not merged into their
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

Any other issues, file a [bug report on
github](https://github.com/kdmurray91/seqhax/issues).

Documentation
=============

The `seqhax` command has many subcommands. The commands, along with a synopsis
of their actions, are displayed when one types `seqhax` with no arguments. 
At the time of writing, these were:

```
$ seqhax
USAGE:
    seqhax PROGRAM [options]

where PROGRAM is one of:
    randfq  -- Generate a random sequence file
    filter  -- Filter reads from a sequence file
    seq     -- Miscellaneous sequence modification
```

The usage of each subcommand can be obtained using the `-h` flag to that
parameter, e.g. `seqhax seq -h`.


## Sub-commands

#### `randfq`

Generates a fasta or fastq file containing sequences with random sequences.

#### `filter`

Removes sequences based on certain criteria:

- Length
* Pairing


#### `seq`

Implements the following functions

- Addition of a constant prefix and/or suffix to each sequence.


License
=======

GPL v3 (see ./LICENSE).

Copyright (c) 2014-2016 Kevin Murray.
