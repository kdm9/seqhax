#!/usr/bin/env python3
import pytest
import subprocess as sp
import os
import shutil
import shlex
import hashlib


def setup(config=None):
    pass


def teardown(config=None):
    sp.check_call("rm -f _*", shell=True)


def run(*args):
    if len(args) == 1:
        args = shlex.split(args[0])
    cmd = ["seqhax", "pairs", *args]
    sp.check_call(cmd, shell=False)


def sha(path):
    h = hashlib.sha512()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def check_eq(f1, f2):
    assert(sha(f1) == sha(f2))


def check_empty(path):
    with open(path, "rb") as fh:
        assert(len(fh.read()) == 0)


def test_pairs():
    run("-f -1 _r1 -2 _r2 -u _rs data/paired_il.fq")
    check_eq("_r1", "data/paired_r1.fq")
    check_eq("_r2", "data/paired_r2.fq")
    check_empty("_rs")
