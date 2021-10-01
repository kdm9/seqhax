#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <regex>
#include "kmseq.hh"


using namespace std;
using namespace kmseq;

struct RebarcodeOptions {
    string r1file;
    string r2file;
    string outfile;
};

bool extract_name_idxseq(const string &header, string &name, vector<string> &idxseqs) {
    regex re("(\\S+) [^:]+:[^:]+:[^:]+:([^+:]+)(\\+([^:+]+))?");
    smatch m;

    idxseqs.clear();

    if (regex_match(header, m, re)) {
        name = m[1].str();
        idxseqs.emplace_back(m[2].str());
        if (m.size() == 5) {
            idxseqs.emplace_back(m[4].str());
        } else {
            idxseqs.emplace_back();
        }
        return true;
    }
    name = "";
    return false;
}


void output_rp(ostream &out, const KSeqPair &sp, const string &id1, const string &id2)
{
    out << "@" << sp.r1.name << endl;
    out << id1 << sp.r1.seq << endl;
    out << "+" << endl;
    out << string(id1.size(), 'I') << sp.r1.qual << endl;

    out << "@" << sp.r2.name << endl;
    out << id2 << sp.r2.seq << endl;
    out << "+" << endl;
    out << string(id2.size(), 'I') << sp.r2.qual << endl;
}

inline void helpmsg(void)
{
    cerr << "USAGE:" << endl;
    cerr << "    seqhax rebarcode [-o OUTPUT] r1 r2" << endl
         << endl;
    cerr << "OPTIONS:"<< endl;
    cerr << "    -o FILE    Output interleaved reads to FILE. Use - for stdout. (default: no output)"<< endl;
}

int parse_args(RebarcodeOptions &opt, int argc, char *argv[])
{

    int c;
    int ret = EXIT_SUCCESS;
    opt.outfile = "-";
    while ((c = getopt(argc, argv, "o:")) > 0) {
        switch (c) {
            case 'o':
                opt.outfile = optarg;
                break;
            case '?':
                ret = EXIT_FAILURE;
            case 'h':
                helpmsg();
                return ret;
        }
    }

    if (optind + 1 >= argc) {
        helpmsg();
        return EXIT_FAILURE;
    }

    opt.r1file = argv[optind++];
    opt.r2file = argv[optind++];
    if (opt.outfile == "-") opt.outfile = "/dev/fd/1";
    return 0;
}

int rebarcode_main(int argc, char *argv[])
{
    RebarcodeOptions opt;

    if (parse_args(opt, argc, argv) != 0) {
        return 1;
    }
    ofstream output;
    if (opt.outfile == "-") {
        output.open("/dev/stdout");
    } else {
        output.open(opt.outfile);
    }

    cerr << "From " << opt.r1file << " and " << opt.r2file << " to " << opt.outfile << endl;

    KSeqPairReader seqs(opt.r1file, opt.r2file);
    size_t count = 0;
    for (KSeqPair sp; seqs.next_pair(sp); count++) {
        string n1, n2;
        vector<string> id1, id2;
        extract_name_idxseq(sp.r1.name, n1, id1);
        extract_name_idxseq(sp.r2.name, n2, id2);
        if (n1 != n2) {
            cerr << "Pair mismatch! Are you sure these are paired-end reads?" << endl;
            return 1;
        }
        if (id1[0] != id2[0] || id1[1] != id2[1]) {
            output_rp(output, sp, "", "");
        } else {
            output_rp(output, sp, id1[0], id1[1]);
        }
    }
    cerr << "Processed " << count << " pairs" << endl;
    return 0;
}

#ifdef REBARCODE_STANDALONE
int main(int argc, char *argv[])
{
    return rebarcode_main(argc, argv);
}
#endif
