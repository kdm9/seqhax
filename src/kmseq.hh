// A thin C++ wrapper around kseq.h
//
// Copyright (c) 2017 Kevin Murray <kdmfoss@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef KMSEQ_HH_T0ZK1CKV
#define KMSEQ_HH_T0ZK1CKV

#include <zlib.h>
#include <string>
#include <ostream>
#include <vector>
#include <stdexcept>
#include "kseq.h"

namespace kmseq
{

using namespace std;

KSEQ_INIT(gzFile, gzread)

struct KSeq
{
    string name;
    string seq;
    string qual;
};


class KSeqReader
{
public:
    KSeqReader()
        : _fp(nullptr)
        , _seq(nullptr)
    {
    }
    KSeqReader(const string &filename)
    {
        this->open(filename);
    }

    ~KSeqReader()
    {
        kseq_destroy(_seq);
        gzclose(_fp);
    }

    void open(const string &filename)
    {
        _fp = gzopen(filename.c_str(), "r");
        if (_fp == NULL) {
            throw runtime_error(string("Could not open file: ") + filename);
        }
        _seq = kseq_init(_fp);
    }

    bool next_read(KSeq &ks)
    {
        if (_seq == nullptr) return false;
        int l = kseq_read(_seq);
        if (l < 1) return false;
        ks.name.assign(_seq->name.s, _seq->name.l);
        ks.seq.assign(_seq->seq.s, _seq->seq.l);
        ks.qual.assign(_seq->qual.s, _seq->qual.l);
        return true;
    }

    size_t next_chunk(vector<KSeq> &sequences, size_t max)
    {
        size_t count = 0;
        for (; count < max; count++) {
            if (count >= sequences.size()) {
                sequences.emplace_back();
            }
            if (!this->next_read(sequences.back())) break;
        }
        sequences.resize(count);
        return count;
    }

protected:
    gzFile _fp;
    kseq_t *_seq;
};

struct KSeqPair
{
    KSeq r1;
    KSeq r2;
};

class KSeqPairReader
{
public:
    KSeqPairReader()
        : _r1()
        , _r2()
        , _interleaved(false)
    {
    }
    KSeqPairReader(const string &r1file, const string &r2file)
    {
        this->open(r1file, r2file);
    }

    KSeqPairReader(const string &interleavedfile)
    {
        this->open(interleavedfile);
    }

    void open(const string &r1file, const string &r2file)
    {
        _r1.open(r1file);
        _r2.open(r2file);
        _interleaved = false;
    }

    void open(const string &interleavedfile)
    {
        _r1.open(interleavedfile);
        _interleaved = true;
    }

    bool next_pair(KSeqPair &ks)
    {
        if (_interleaved) {
            return _r1.next_read(ks.r1) && _r1.next_read(ks.r2);
        } else {
            return _r1.next_read(ks.r1) && _r2.next_read(ks.r2);
        }
    }

    size_t next_chunk(vector<KSeqPair> &pairs, size_t max)
    {
        size_t count = 0;
        for (; count < max; count++) {
            if (count >= pairs.size()) {
                pairs.emplace_back();
            }
            if (!this->next_pair(pairs.back())) break;
        }
        pairs.resize(count);
        return count;
    }

protected:
    KSeqReader _r1;
    KSeqReader _r2;
    bool _interleaved;
};


} /* end namespace kmseq */

inline std::ostream &operator<<(std::ostream &out, const kmseq::KSeq &r)
{ 
    using std::endl;
    out << "@" << r.name << endl;
    out << r.seq << endl;
    out << "+" << endl;
    out << r.qual << endl;
    return out;
}

inline std::ostream &operator<<(std::ostream &out, const kmseq::KSeqPair &rp)
{ 
    out << rp.r1;
    out << rp.r2;
    return out;
}

#endif /* end of include guard: KMSEQ_HH_T0ZK1CKV */

// vim:set et sw=4 ts=4:
