#!/usr/bin/env python3
from __future__ import print_function


MAXSIZE = 10000000

seq = ""
for i in range(1, MAXSIZE + 1):
    seq += 'A'
    print("@", i, "\n", seq, "\n+\n", seq, '\n', sep='', end='')
