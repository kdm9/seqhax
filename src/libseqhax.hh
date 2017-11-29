#ifndef LIBSEQHAX_HH_KWLGOZPY
#define LIBSEQHAX_HH_KWLGOZPY

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "kmseq.hh"

namespace seqhax
{

using namespace std;
using namespace kmseq;

void extract_readname(const string &header, string &name)
{
    size_t i = header.find(" ");
    if (i == string::npos) {
        i = header.size();
    }
    // @asldfkjaskdf/1<END_OF_STRING>
    // @asldfkjaskdf/1 1:N:blah
    //                ^
    //                i is here
    if (header[i-2] == '/' && (header[i-1] == '1' || header[i-1] == '2')) {
        i -= 2;
    }
    name = header.substr(0, i);
}


} /* seqhax */

#endif /* end of include guard: LIBSEQHAX_HH_KWLGOZPY */