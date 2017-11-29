#include <iostream>
#include <fstream>
#include <iterator>
#include <string>

#include <boost/program_options.hpp>

#include "kmseq.hh"
#include "libseqhax.hh"

using namespace std;
using namespace kmseq;
using namespace seqhax;

struct PECheckOptions {
    string r1file;
    string r2file;
    string outfile;
};


int parse_args(PECheckOptions &opt, int argc, char *argv[])
{
    namespace po = boost::program_options;
    po::options_description flags("Options");
    flags.add_options()
        ("help,h",
            "Print this help")
        ("outfile,o", po::value<string>(&opt.outfile),
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
    if (opt.outfile == "-") opt.outfile = "/dev/fd/1";
    return 0;
}


int pecheck_main(int argc, char *argv[])
{
    PECheckOptions opt;
    ofstream output;

    if (parse_args(opt, argc, argv) != 0) {
        return 1;
    }
    cerr << "From " << opt.r1file << " and " << opt.r2file << " to " << opt.outfile << endl;

    if (opt.outfile != "") {
        output.open(opt.outfile);
    }

    KSeqPairReader seqs(opt.r1file, opt.r2file);
    size_t count = 0;
    for (KSeqPair sp; seqs.next_pair(sp); count++) {
        string n1, n2;
        extract_readname(sp.r1.name, n1);
        extract_readname(sp.r2.name, n2);
        if (n1 != n2) {
            cerr << "Pair mismatch! Are you sure these are paired-end reads?" << endl;
            cerr << endl;
            cerr << "At pair " << count << endl;
            cerr << "'" << n1 << "' != '" << n2 << "'" << endl;
            return 1;
        }
        if (output.is_open()) {
            output << sp;
        }
    }
    cerr << "Processed " << count << " pairs" << endl;
    return 0;
}

#ifdef PECHECK_STANDALONE
int main(int argc, char *argv[])
{
    return pecheck_main(argc, argv);
}
#endif
