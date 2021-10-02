#include <iostream>

#include <map>
#include <cmath>
#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

using namespace std;

int percent(float val) {
    return round(val * 100.0);
}

void usage() {
    cerr 
        << "USAGE:" << endl
        << "    htshax bcfhist BCF_FILE" << endl
        << endl
        << "Calculates histograms of allele frequency and missingness across SNPs."
        << endl;
}


int bcfhist_main(int argc, char **argv) {
    // TODO: expose sample list and regions to CLI.

    if (argc != 2) {
        usage();
        return EXIT_FAILURE;
    }
    string filename = argv[1];

    if (filename == "--help" || filename == "-h") {
        usage();
        return EXIT_FAILURE;
    }

    string region;
    vector<string> samples;
    string samplecsv;
    bool success = true;
    map<int, uint64_t> af_hist;
    map<int, uint64_t> miss_hist;
    map<int, uint64_t> qual_hist;
    map<int, uint64_t> dp_hist;

    for (int i=0; i<=100; i++) {
        af_hist[i] = 0;
        miss_hist[i] = 0;
    }

    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    if (success && samplecsv.size() > 0) {
        if (bcf_sr_set_samples(sr, samplecsv.c_str(), 0) != 1) {
            cerr << "Failed to set samples for " << filename << ": " << samplecsv << endl;
            success = false;
        }
    }
    if (success && region.size() > 0) {
        if (bcf_sr_set_regions(sr, region.c_str(), 0) != 0) {
            cerr << "Failed to set region for " << filename << endl;
            success = false;
        }
    }
    if (success && bcf_sr_add_reader(sr, filename.c_str()) != 1) {
        cerr << "Failed to add reader for " << filename << endl;
        success = false;
    }
    bcf_hdr_t *header = sr->readers[0].header;

    int32_t *buffer = NULL;
    int32_t buffersz = 0;
    while (success && bcf_sr_next_line(sr)) {
        bcf1_t *record = bcf_sr_get_line(sr,0);
        int32_t nentry = bcf_get_format_int32(header, record, "GT", &buffer, &buffersz);
        int nmiss = 0;
        int nalt = 0;
        if (nentry < 0) {
            success = false;
        }
        for (int32_t i=0; i<nentry; i += 2) {
            int32_t a=bcf_gt_allele(buffer[i]), b=bcf_gt_allele(buffer[i+1]);
            if (a < 0) nmiss++; else nalt += a;
            if (b < 0) nmiss++; else nalt += b;
        }
        float miss = (float)nmiss / nentry;
        float af = (float)nalt / (nentry - nmiss);
        bcf_info_t *dp_info = bcf_get_info(header, record, "DP");
        float dp = -1;
        if (dp_info != NULL) {
            assert(dp_info->len == 1);
            switch (dp_info->type) {
                case BCF_BT_INT8:
                case BCF_BT_INT16:
                case BCF_BT_INT32:
                case BCF_BT_INT64:
                    dp = dp_info->v1.i;
                    break;
                case BCF_BT_FLOAT:
                    dp = dp_info->v1.f;
                    break;
            }
        }


        //cerr << bcf_hdr_id2name(header, record->rid) << "\t" << record->pos << "\t" << nmiss << "\t" << nalt << "\t" << nentry - nmiss << "\t" << af << "\t" << miss << endl;

        //af_hist(af);
        //mis_hist(miss);
        af_hist[percent(af)]++;
        miss_hist[percent(miss)]++;
        qual_hist[(int)round(record->qual)]++;
        dp_hist[(int)round(dp)]++;
    }
    if (sr->errnum) {
        cerr << "Error: " << bcf_sr_strerror(sr->errnum) << endl;
    }

    if (success) {
        cout << "metric\tpercent\tnsnp" << endl;
        for (const auto &[key, value]: af_hist) {
            cout << "af\t" << key << "\t" << value  << endl;
        }
        for (const auto &[key, value]: miss_hist) {
            cout << "miss\t" << key << "\t" << value  << endl;
        }
        for (const auto &[key, value]: qual_hist) {
            cout << "qual\t" << key << "\t" << value  << endl;
        }
        for (const auto &[key, value]: dp_hist) {
            cout << "dp\t" << key << "\t" << value  << endl;
        }
    }

    if (buffer != NULL) {
        free(buffer);
        buffer = NULL;
        buffersz = 0;
    }

    bcf_sr_destroy(sr);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

