#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <regex>
#include "kmseq.hh"

#include <boost/program_options.hpp>

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


int parse_args(RebarcodeOptions &opt, int argc, char *argv[])
{
    namespace po = boost::program_options;
    po::options_description flags("Options");
    flags.add_options()
        ("help,h",
            "Print this help")
        ("outfile,o", po::value<string>(&opt.outfile)->default_value("/dev/stdout"),
            "Output filename")
        ("r1", po::value<string>(&opt.r1file)->required(),
            "First read filename")
        ("r2", po::value<string>(&opt.r2file)->required(),
            "Second read filename");
    po::positional_options_description pos;
    pos.add("r1", 1);
    pos.add("r2", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
              .options(flags).positional(pos).run(), vm);

    if (vm.count("help")) {
        cout << flags << endl;
        exit(0);
    }
    try {
        po::notify(vm);
    } catch (po::error &e) {
        cerr << e.what() << endl << endl;
        cerr << flags << endl;
        return 1;
    }
    return 0;
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
