#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <unistd.h>


int clihist_main(int argc, char **argv)
{
    using namespace std;
    // Makes cin not horribly slow. See https://stackoverflow.com/questions/9371238
    cin.sync_with_stdio(false);

    if (isatty(fileno(stdout))) {
        cerr << "clihist expects input on stdin and takes no args. Just pipe into me!" << endl << endl;
        cerr << "WARNING: clihist should not be used with more than a few thousand possible values (beware floating point numbers)!!" << endl;
    }

    map<string, size_t> histo;
    string line;
    while (getline(cin, line)) {
        if (cin.eof()) break;
        ++histo[line];
    };

    for (const auto &h: histo) {
        cout << h.first << "\t" << h.second << "\n";
    }
    return EXIT_SUCCESS;
}
