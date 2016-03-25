/*
 * About the implemantation of a HMM using a nonparametric distribution
 * to model the hidden variable P (i.e., underlying methylation ratio),
 * conditioned on which is the binomial distribution with a fixed N for
 * each CpG site.
 *
 * Author:  Shiqi Tu
 * Version: 1.0
 * Date:    03-09-2016
 * 
 */

#ifndef NONPARAMETRIC_RATIO_HMM_HPP
#define NONPARAMETRIC_RATIO_HMM_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <iostream>

class CpG_block {

public:
    CpG_block(): start(0), end(0), methy_c(0), coverage(0) {}
    CpG_block(std::string c, size_t pos, size_t i, size_t j):
              chr(c), start(pos), end(pos + 1), methy_c(i), coverage(j) {}
    CpG_block(std::string c, size_t s, size_t e, size_t i, size_t j):
              chr(c), start(s), end(e), methy_c(i), coverage(j) {}
    void set(std::string c, size_t s, size_t e, size_t i, size_t j) {
        chr = c; start = s; end = e; methy_c = i; coverage = j;
    }
    void set(std::string c, size_t pos, size_t i, size_t j) {
        set(c, pos, pos + 1, i, j);
    }
    
    std::string chr;
    size_t start, end;
    size_t methy_c;
    size_t coverage;
};

typedef std::vector<CpG_block> obs_seq;
typedef std::vector<double> log_prob_dist;
typedef std::vector<double> prob_dist;

class NonparametricRatioHMM {

public:
    typedef double log_likelihood;
    NonparametricRatioHMM(size_t ns, size_t nb=40,
                          double mp=0.0, double tol=8.0E-5, size_t mi=60);
    
    void setParas(const prob_dist &pi, const std::vector<prob_dist> &trans,
                  const std::vector<prob_dist> &emission);
    void getParas(prob_dist &pi, std::vector<prob_dist> &trans,
                  std::vector<prob_dist> &emission) const;
    
    // train parameters & return the converged log-likelihood for the observations.
    log_likelihood trainParas(const std::vector<obs_seq> &obs, const prob_dist &init_pi,
                              const std::vector<prob_dist> &init_trans,
                              const std::vector<prob_dist> &init_emission,
                              bool verbose=false, std::ostream &os=std::cout);
    
    log_likelihood getLogLikelihood(const std::vector<obs_seq> &obs) const;
    static bool is_log0(log_likelihood x);
    
    // posterior decoding for hidden states.
    log_likelihood posterior_prob(const std::vector<obs_seq> &obs,
                                  std::vector<std::vector<prob_dist> > &post_dist_states) const;

private:
    static double log_sum_exp(const std::vector<double> &x, double cutoff);
    static double log_sum_exp(double x, double y, double cutoff);
    
    void adjustParas(log_prob_dist &x, double cutoff, double mark);
    static void switchMarksUtil(log_prob_dist &x, double cutoff, double mark);
    void switchMarks(double cutoff, double mark);
    
    void p2logp(const prob_dist &x, log_prob_dist &y);
    static void logp2p(const log_prob_dist &x, prob_dist &y);
    
    static void log_choose(const std::vector<obs_seq> &obs, std::vector<std::vector<double> > &ret);
    typedef std::vector<std::vector<double> > Matrix;
    void get_site_log_emission(std::vector<Matrix> &ret, const std::vector<obs_seq> &obs,
                               const std::vector<std::vector<double> > &log_binomial_coeff,
                               const std::vector<double> &log_ratio, const double cutoff) const;
    void forward(std::vector<Matrix> &alpha, const std::vector<Matrix> &site_log_emission,
                 const double cutoff) const;
    void backward(std::vector<Matrix> &gama, std::vector<log_prob_dist> &n_log_trans,
                  log_prob_dist &n_log_pi, const std::vector<Matrix> &alpha,
                  const std::vector<Matrix> &site_log_emission, const double cutoff) const;
    void fit_emission(std::vector<log_prob_dist> &n_log_emission, const std::vector<Matrix> &gama,
                      const std::vector<Matrix> &site_log_emission, const std::vector<obs_seq> &obs,
                      const std::vector<std::vector<double> > &log_binomial_coeff,
                      const std::vector<double> &log_ratio, const double cutoff) const;
    
    static log_likelihood getLogLikelihoodUtil(const std::vector<Matrix> &alpha,
                                               const double cutoff);
    void printInfo(size_t iter_n, log_likelihood logP, std::ostream &os,
                   const double cutoff) const;
    
    // log of minimum allowed value for the "parameter" probabilities that are not 0.
    double log_min_prob;
    
    // for training parameters.
    double log_likelihood_tol;
    size_t max_iterations;
    
    // state == 0 means parameters haven't been set or trained yet.
    int state;
    
    // parameters defining the HMM.
    const size_t N, nbins;
    log_prob_dist log_pi;
    std::vector<log_prob_dist> log_trans;
    std::vector<log_prob_dist> log_emission;
    
    // special mark: for a probability p, p == 0 <=> log(p) >= MARK
    static const double MARK;
};

#endif

