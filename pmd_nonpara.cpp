/*
 * This program is for recognizing partially methylated domains (PMDs)
 * in a given methylome using a nonparametric-ratio HMM.
 *
 * Author:  Shiqi Tu
 * Version: 1.0
 * Date:    03-24-2016
 *
 */


#include "NonparametricRatioHMM.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <map>

using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::runtime_error;
using std::pair;

typedef vector<pair<long, long> > block_seq;

class CpG_site {

public:
    CpG_site(const string &line);
    void set(const string &line);
    
    string chr;
    size_t pos;
    size_t methy_c, coverage;
};

CpG_site::CpG_site(const string &line) {
    set(line);
}

void CpG_site::set(const string &line) {
    std::istringstream is(line);
    string dump;
    double ratio = 0.0;
    is >> chr >> pos >> dump >> dump >> ratio >> coverage;
    if (!is) {
        throw runtime_error("Invalid line format for describing a cytosine site!");
    }
    if (ratio < 0.0 || ratio > 1.0) {
        throw runtime_error("Methylation level must be in the range of [0, 1] !");
    }
    methy_c = round(coverage * ratio);
}

// perform a segmentation of the methylome into multiple observation sequences.
void segment_methylome(std::istream &meth, size_t bin_size, size_t max_gap,
                       vector<obs_seq> &obs, vector<block_seq> &to_trim) {
    obs.clear();
    to_trim.clear();
    
    string line;
    if (!std::getline(meth, line)) return;
    
    CpG_site cpg(line);
    CpG_block curr_bin(cpg.chr, cpg.pos, cpg.pos + bin_size, cpg.methy_c, cpg.coverage);
    size_t pre_cpg = cpg.pos;
    pair<long, long> fl(cpg.pos, -1);
    
    obs.push_back(obs_seq());
    vector<obs_seq>::iterator curr_seq = obs.end() - 1;
    to_trim.push_back(block_seq());
    vector<block_seq>::iterator curr_trim = to_trim.end() - 1;
    
    while (std::getline(meth, line)) {
        cpg.set(line);
        if (cpg.chr != curr_bin.chr || cpg.pos > pre_cpg + max_gap) {
            curr_seq -> push_back(curr_bin);
            fl.second = pre_cpg; curr_trim -> push_back(fl);
            curr_bin.set(cpg.chr, cpg.pos, cpg.pos + bin_size, 0, 0);
            fl.first = cpg.pos;
            obs.push_back(obs_seq());
            curr_seq = obs.end() - 1;
            to_trim.push_back(block_seq());
            curr_trim = to_trim.end() - 1;
        }
        else if (cpg.pos >= curr_bin.end) {
            curr_seq -> push_back(curr_bin);
            fl.second = pre_cpg; curr_trim -> push_back(fl);
            curr_bin.start = curr_bin.end; curr_bin.end += bin_size;
            curr_bin.methy_c = curr_bin.coverage = 0;
            fl.first = fl.second = -1;
            while (cpg.pos >= curr_bin.end) {
                curr_seq -> push_back(curr_bin);
                curr_trim -> push_back(fl);
                curr_bin.start = curr_bin.end; curr_bin.end += bin_size;
            }
            fl.first = cpg.pos;
        }
        pre_cpg = cpg.pos;
        curr_bin.methy_c += cpg.methy_c;
        curr_bin.coverage += cpg.coverage;
    }
    curr_seq -> push_back(curr_bin);
    fl.second = pre_cpg; curr_trim -> push_back(fl);
}

// perform an empirical estimation of the emission parameters for
// the initialization of HMM.
void initialize(const vector<obs_seq> &obs, prob_dist &pi, vector<prob_dist> &trans,
                vector<prob_dist> &emi, size_t nbins, size_t bin_size) {
    // state 1 refers to being in a PMD.
    pi.resize(2, 1.0);
    
    trans.resize(2);
    trans[0].resize(2, 1.0);
    trans[1].resize(2, 1.0);
    
    emi.resize(2);
    emi[0].resize(nbins + 1, 1.0);
    emi[1].resize(nbins + 1, 1.0);
    
    size_t range = 20000, min_cover = 10;
    double pmd_low = 0.15, pmd_up = 0.6;
    double background_low = 0.75;
    double min_sigma = 0.03;
    
    size_t num = range / bin_size;
    if (range % bin_size != 0) ++num;
    obs_seq con_bins;
    for (vector<obs_seq>::const_iterator obs_iter = obs.begin();
         obs_iter != obs.end(); ++obs_iter) {
        if (obs_iter -> size() < num) continue;
        size_t i = 0, total_methy_c = 0, total_coverage = 0;
        con_bins.clear();
        while (i < num) {
            con_bins.push_back(obs_iter -> at(i));
            total_methy_c += con_bins[i].methy_c;
            total_coverage += con_bins[i].coverage;
            ++i;
        }
        size_t front = num - 1, tail = 0;
        while (true) {
            double ratio = (total_coverage == 0) ? 0.0 : static_cast<double>(total_methy_c) / total_coverage;
            if (total_coverage < min_cover || ratio < pmd_low || (ratio > pmd_up && ratio < background_low)) {
                if (i == obs_iter -> size()) break;
                total_methy_c -= con_bins[tail].methy_c;
                total_coverage -= con_bins[tail].coverage;
                tail = (tail + 1) % num; front = (front + 1) % num;
                con_bins[front] = obs_iter -> at(i++);
                total_methy_c += con_bins[front].methy_c;
                total_coverage += con_bins[front].coverage;
                continue;
            }
            int flag = (ratio >= background_low) ? 0 : 1;
            for (size_t j = 0; j != num; ++j) {
                if (con_bins[j].coverage < min_cover) continue;
                util_emission_update(emi[flag], con_bins[j].methy_c, con_bins[j].coverage, min_sigma);
            }
            
        }
        
    }
    
}















// for now,
// pmd_nonpara -o fitting_process input_fn
int main(int argc, char **argv) {
    std::ifstream meth(argv[3]);
    if (!meth) {
        cerr << "Can't open the input file '" << argv[3] << "'!" << endl;
        return -1;
    }
    size_t bin_size = 500, max_gap = 1000;
    
    // segmentation.
    vector<obs_seq> obs;
    vector<block_seq> to_trim;
    segment_methylome(meth, bin_size, max_gap, obs, to_trim);
    meth.close();
    
    // fitting a nonparametric-ratio HMM.
    std::ofstream outf(argv[2]);
    if (!outf) {
        cerr << "Can't create the output file '" << argv[2] << "'!" << endl;
        return -2;
    }
    
    // initialize the parameters.
    const size_t nbins = 40;
    prob_dist pi;
    vector<prob_dist> trans, emi;
    initialize(obs, pi, trans, emi, nbins, bin_size);
    
    NonparametricRatioHMM hmm(2, nbins);
    hmm.trainParas(obs, pi, trans, emi, true, outf);
    outf.close();
    
    return 0;
}

