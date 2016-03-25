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

#include "NonparametricRatioHMM.hpp"
#include <cmath>
#include <stdexcept>

using std::vector;
using std::runtime_error;
using std::endl;

const double NonparametricRatioHMM::MARK = 1.0;

/////////////////////////////////////////////////////////////////////////////////////////////////
// For the following functions, log(p) >= cutoff means p == 0 and vice versa.
// In fact, log(p) == 2 * cutoff when p == 0, and log(p) <= cutoff / 2 otherwise.
// cutoff is supposed to be greater than MARK.
// For the final "parameter" probabilities, log(p) == 2 * MARK when p == 0.
// In this case, set mark to be 2 * MARK.

NonparametricRatioHMM::NonparametricRatioHMM(size_t ns, size_t nb, double mp,
                                             double tol, size_t mi): log_likelihood_tol(tol),
                                             max_iterations(mi), state(0), N(ns), nbins(nb),
                                             log_pi(N, 0.0), log_trans(N, log_pi),
                                             log_emission(N, log_prob_dist(nbins + 1, 0.0)) {
    if (mp > 1.0) {
        throw runtime_error("Minimum allowed probability for each parameter is > 1 !");
    }
    else if (mp <= 0.0) {
        log_min_prob = MARK * 2.0;
    }
    else {
        log_min_prob = log(mp);
    }
    
    if (N == 0) {
        throw runtime_error("There should be at least 1 underlying state!");
    }
    if (nbins == 0) {
        throw runtime_error("There should be at least 1 bin to segment the probability range [0, 1] !");
    }
}

double NonparametricRatioHMM::log_sum_exp(const std::vector<double> &x, double cutoff) {
    double lse = cutoff * 2.0, temp = 0.0;
    size_t i = 0;
    for (; i != x.size(); ++i) {
        if (x[i] < cutoff) {
            lse = x[i++];
            break;
        }
    }
    for (; i != x.size(); ++i) {
        temp = x[i];
        if (temp < cutoff) {
            if (lse < temp) {
                lse = temp + log(1.0 + exp(lse - temp));
            }
            else {
                lse = lse + log(1.0 + exp(temp - lse));
            }
        }
    }
    return lse;
}

double NonparametricRatioHMM::log_sum_exp(double x, double y, double cutoff) {
    if (x >= cutoff) return y;
    if (y >= cutoff) return x;
    if (x < y) {
        double temp = x;
        x = y;
        y = temp;
    }
    return x + log(1.0 + exp(y - x));
}

void NonparametricRatioHMM::adjustParas(log_prob_dist &x, double cutoff, double mark) {
    double temp = cutoff + 2.0 * log(x.size());
    switchMarksUtil(x, cutoff, temp * 2.0);
    cutoff = temp;
    
    double lse = log_sum_exp(x, cutoff);
    if (lse >= cutoff) {
        throw runtime_error("None of the probabilities in the given distribution is greater than 0!");
    }
    bool flag = false;
    for (size_t i = 0; i != x.size(); ++i) {
        if (x[i] >= cutoff) {
            x[i] = mark;
        }
        else {
            x[i] -= lse;
            if (log_min_prob < MARK && x[i] < log_min_prob) {
                x[i] = mark;
                flag = true;
            }
        }
    }
    if (!flag) return;
    
    cutoff = mark / 2.0;
    lse = log_sum_exp(x, cutoff);
    if (lse >= cutoff) {
        throw runtime_error("None of the probabilities in the given distribution is greater than 0!");
    }
    for (size_t i = 0; i != x.size(); ++i) {
        if (x[i] < cutoff) x[i] -= lse;
    }
}

void NonparametricRatioHMM::switchMarksUtil(log_prob_dist &x, double cutoff, double mark) {
    for (size_t i = 0; i != x.size(); ++i) {
        if (x[i] >= cutoff) {
            x[i] = mark;
        }
    }
}

void NonparametricRatioHMM::switchMarks(double cutoff, double mark) {
    switchMarksUtil(log_pi, cutoff, mark);
    for (size_t i = 0; i != N; ++i) {
        switchMarksUtil(log_trans[i], cutoff, mark);
        switchMarksUtil(log_emission[i], cutoff, mark);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void NonparametricRatioHMM::p2logp(const prob_dist &x, log_prob_dist &y) {
    double p = 0.0, max = MARK / 2.0;
    for (size_t i = 0; i != y.size(); ++i) {
        p = x.at(i);
        if (p < 0.0) {
            throw runtime_error("Probabilities can't be negative!");
        }
        else if (p > 0.0) {
            p = log(p);
            max = p > max ? p : max;
        }
    }
    double cutoff = max * 2.0;
    for (size_t i = 0; i != y.size(); ++i) {
        p = x.at(i);
        if (p > 0.0) {
            y[i] = log(p);
        }
        else {
            y[i] = cutoff * 2.0;
        }
    }
    adjustParas(y, cutoff, MARK * 2.0);
}

void NonparametricRatioHMM::logp2p(const log_prob_dist &x, prob_dist &y) {
    y.resize(x.size());
    for (size_t i = 0; i != x.size(); ++i) {
        y[i] = x[i] < MARK ? exp(x[i]) : 0.0;
    }
}

void NonparametricRatioHMM::setParas(const prob_dist &pi, const std::vector<prob_dist> &trans,
                                     const std::vector<prob_dist> &emission) {
    p2logp(pi, log_pi);
    for (size_t i = 0; i != N; ++i) {
        p2logp(trans.at(i), log_trans.at(i));
        p2logp(emission.at(i), log_emission.at(i));
    }
    state = 1;
}

void NonparametricRatioHMM::getParas(prob_dist &pi, std::vector<prob_dist> &trans,
                                     std::vector<prob_dist> &emission) const {
    if (state == 0) {
        throw runtime_error("HMM parameters haven't been set or trained yet!");
    }
    
    logp2p(log_pi, pi);
    trans.resize(N);
    emission.resize(N);
    for (size_t i = 0; i != N; ++i) {
        logp2p(log_trans.at(i), trans.at(i));
        logp2p(log_emission.at(i), emission.at(i));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void NonparametricRatioHMM::log_choose(const std::vector<obs_seq> &obs,
                                       std::vector<std::vector<double> > &ret) {
    ret.resize(obs.size());
    vector<double>::iterator iter;
    for (size_t i = 0; i != obs.size(); ++i) {
        ret[i].resize(obs[i].size());
        iter = ret[i].begin();
        for (obs_seq::const_iterator iter1 = obs[i].begin();
             iter1 != obs[i].end(); ++iter1) {
            size_t a = iter1 -> methy_c, b = iter1 -> coverage;
            if (a > b) {
                throw runtime_error("Number of methylated reads can't exceed the corresponding coverage!");
            }
            double res = 0.0;
            while (a != 0) {
                res += log(b--) - log(a--);
            }
            *iter++ = res;
        }
    }
}

void NonparametricRatioHMM::get_site_log_emission(std::vector<Matrix> &ret, const std::vector<obs_seq> &obs,
                                                  const std::vector<std::vector<double> > &log_binomial_coeff,
                                                  const std::vector<double> &log_ratio, const double cutoff) const {
    vector<double> temp(nbins + 1);
    for (size_t i = 0; i != obs.size(); ++i) {
        for (size_t j = 0; j != obs[i].size(); ++j) {
            vector<double> *site = &ret[i][j];
            size_t a = obs[i][j].methy_c, b = obs[i][j].coverage - a;
            double log_coeff = log_binomial_coeff[i][j];
            for (size_t k = 0; k != N; ++k) {
                for (size_t bin = 0; bin <= nbins; ++bin) {
                    temp[bin] = 2.0 * cutoff;
                    double logP = log_emission[k][bin];
                    if (logP >= cutoff) continue;
                    double lp1 = log_ratio[bin], lp2 = log_ratio[nbins - bin];
                    if (a != 0 && lp1 >= cutoff) continue;
                    lp1 *= a;
                    if (b != 0 && lp2 >= cutoff) continue;
                    lp2 *= b;
                    temp[bin] = log_coeff + lp1 + lp2 + logP;
                }
                site -> at(k) = log_sum_exp(temp, cutoff);
            }
        }
    }
}

void NonparametricRatioHMM::forward(std::vector<Matrix> &alpha, const std::vector<Matrix> &site_log_emission,
                                    const double cutoff) const {
    vector<double> temp(N);
    for (size_t i = 0; i != alpha.size(); ++i) {
        Matrix::iterator pre = alpha[i].begin();
        Matrix::const_iterator emi = site_log_emission[i].begin();
        for (size_t j = 0; j != N; ++j) {
            double lp1 = log_pi[j], lp2 = emi -> at(j);
            pre -> at(j) = (lp1 >= cutoff || lp2 >= cutoff) ? (2.0 * cutoff) : (lp1 + lp2);
        }
        ++emi;
        for (Matrix::iterator nxt = pre + 1; nxt != alpha[i].end(); ++nxt) {
            for (size_t j = 0; j != N; ++j) {
                nxt -> at(j) = 2.0 * cutoff;
                double logP = emi -> at(j);
                if (logP >= cutoff) continue;
                for (size_t k = 0; k != N; ++k) {
                    double lp1 = pre -> at(k), lp2 = log_trans[k][j];
                    temp[k] = (lp1 >= cutoff || lp2 >= cutoff) ? (2.0 * cutoff) : (lp1 + lp2);
                }
                double lse = log_sum_exp(temp, cutoff);
                if (lse >= cutoff) continue;
                nxt -> at(j) = logP + lse;
            }
            ++pre;
            ++emi;
        }
    }
}

void NonparametricRatioHMM::backward(std::vector<Matrix> &gama, std::vector<log_prob_dist> &n_log_trans,
                                     log_prob_dist &n_log_pi, const std::vector<Matrix> &alpha,
                                     const std::vector<Matrix> &site_log_emission, const double cutoff) const {
    vector<double> temp(N);
    for (size_t i = 0; i != N; ++i) {
        for (log_prob_dist::iterator iter = n_log_trans[i].begin();
             iter != n_log_trans[i].end(); ++iter) *iter = 2.0 * cutoff;
        n_log_pi[i] = 2.0 * cutoff;
    }
    
    vector<double> *beta_t = new vector<double>(N);
    vector<double> *beta_t_1 = new vector<double>(N);
    for (size_t i = 0; i != alpha.size(); ++i) {
        beta_t -> resize(N, 0.0);
        Matrix::const_reverse_iterator alpha_iter = alpha[i].rbegin(),
                                       emi_iter = site_log_emission[i].rbegin();
        Matrix::reverse_iterator gama_iter = gama[i].rbegin();
        double offset = log_sum_exp(*alpha_iter, cutoff);
        for (size_t j = 0; j != N; ++j) {
            double lp = alpha_iter -> at(j);
            gama_iter -> at(j) = (lp >= cutoff) ? (2.0 * cutoff) : (lp - offset);
        }
        while (++gama_iter != gama[i].rend()) {
            ++alpha_iter;
            for (size_t j = 0; j != N; ++j) {
                for (size_t k = 0; k != N; ++k) {
                    temp[k] = 2.0 * cutoff;
                    double lp = log_trans[j][k];
                    if (lp >= cutoff) continue;
                    double lp1 = emi_iter -> at(k), lp2 = beta_t -> at(k);
                    if (lp1 >= cutoff || lp2 >= cutoff) continue;
                    temp[k] = lp + lp1 + lp2;
                }
                beta_t_1 -> at(j) = log_sum_exp(temp, cutoff);
                double lp1 = alpha_iter -> at(j), lp2 = beta_t_1 -> at(j);
                gama_iter -> at(j) = (lp1 >= cutoff || lp2 >= cutoff) ? (2.0 * cutoff) : (lp1 + lp2 - offset);
                
                // update "n_log_trans".
                if (lp1 >= cutoff) continue;
                for (size_t k = 0; k != N; ++k) {
                    if (temp[k] >= cutoff) continue;
                    n_log_trans[j][k] = log_sum_exp(n_log_trans[j][k], lp1 + temp[k] - offset, cutoff);
                }
            }
            ++emi_iter;
            vector<double> *tmp = beta_t;
            beta_t = beta_t_1;
            beta_t_1 = tmp;
        }
        
        // update "n_log_pi".
        --gama_iter;
        for (size_t j = 0; j != N; ++j) n_log_pi[j] = log_sum_exp(n_log_pi[j], gama_iter -> at(j), cutoff);
    }
    delete beta_t;
    delete beta_t_1;
}

void NonparametricRatioHMM::fit_emission(std::vector<log_prob_dist> &n_log_emission, const std::vector<Matrix> &gama,
                                         const std::vector<Matrix> &site_log_emission, const std::vector<obs_seq> &obs,
                                         const std::vector<std::vector<double> > &log_binomial_coeff,
                                         const std::vector<double> &log_ratio, const double cutoff) const {
    vector<double> temp(nbins + 1);
    for (size_t i = 0; i != N; ++i) {
        for (log_prob_dist::iterator iter = n_log_emission[i].begin();
             iter != n_log_emission[i].end(); ++iter) *iter = 2.0 * cutoff;
    }
    
    for (size_t i = 0; i != obs.size(); ++i) {
        Matrix::const_iterator gama_iter = gama[i].begin(),
                               emi_iter = site_log_emission[i].begin();
        obs_seq::const_iterator obs_iter = obs[i].begin();
        vector<double>::const_iterator lbc_iter = log_binomial_coeff[i].begin();
        while (obs_iter != obs[i].end()) {
            size_t a = obs_iter -> methy_c, b = obs_iter -> coverage - a;
            for (size_t k = 0; k != N; ++k) {
                if (gama_iter -> at(k) >= cutoff) continue;
                for (size_t bin = 0; bin <= nbins; ++bin) {
                    temp[bin] = 2.0 * cutoff;
                    double logP = log_emission[k][bin];
                    if (logP >= cutoff) continue;
                    double lp1 = log_ratio[bin], lp2 = log_ratio[nbins - bin];
                    if (a != 0 && lp1 >= cutoff) continue;
                    lp1 *= a;
                    if (b != 0 && lp2 >= cutoff) continue;
                    lp2 *= b;
                    temp[bin] = (*lbc_iter) + lp1 + lp2 + logP;
                }
                
                // update "n_log_emission".
                for (size_t bin = 0; bin <= nbins; ++bin) {
                    if (temp[bin] >= cutoff) continue;
                    n_log_emission[k][bin] = log_sum_exp(n_log_emission[k][bin],
                                                         gama_iter -> at(k) + temp[bin] - emi_iter -> at(k), cutoff);
                }
            }
            ++gama_iter;
            ++emi_iter;
            ++obs_iter;
            ++lbc_iter;
        }
    }
}

NonparametricRatioHMM::log_likelihood
NonparametricRatioHMM::getLogLikelihoodUtil(const vector<Matrix> &alpha,
                                            const double cutoff) {
    log_likelihood logP = 0.0;
    for (size_t i = 0; i != alpha.size(); ++i) {
        Matrix::const_iterator iter = alpha[i].end() - 1;
        double lse = log_sum_exp(*iter, cutoff);
        if (lse >= cutoff) return 2.0 * cutoff;
        logP += lse;
    }
    return logP;
}

void NonparametricRatioHMM::printInfo(size_t iter_n, log_likelihood logP, std::ostream &os,
                                      const double cutoff) const {
    std::ostream::fmtflags original_flgs = os.flags();
    int prob_pre = 4, logP_pre = ceil(-log(log_likelihood_tol) / log(10.0));
    logP_pre = (logP_pre > 0) ? (logP_pre + 2) : 2;
    int original_pre = os.precision(prob_pre);
    os << std::fixed << std::showpoint;
    
    if (iter_n == 0) os << "Initialization:" << endl;
    else os << "\n================================================================================\n"
            << "After iteration step " << iter_n << ":" << endl;
    
    os << "\nInitial state probability distribution:\npi:";
    for (log_prob_dist::const_iterator iter = log_pi.begin();
         iter != log_pi.end(); ++iter) {
        os << "\t" << ((*iter >= cutoff) ? 0.0 : exp(*iter));
    }
    os << endl;
    
    os << "\nTransition matrix:\nTrans";
    for (size_t i = 0; i != N; ++i) os << "\tstate" << i;
    os << endl;
    for (size_t i = 0; i != N; ++i) {
        os << "state" << i;
        for (log_prob_dist::const_iterator iter = log_trans[i].begin();
             iter != log_trans[i].end(); ++iter) {
            os << "\t" << ((*iter >= cutoff) ? 0.0 : exp(*iter));
        }
        os << endl;
    }
    
    os << "\nEmission matrix:\nRatio";
    for (size_t i = 0; i <= nbins; ++i) os << "\t" << static_cast<double>(i) / nbins;
    os << endl;
    for (size_t i = 0; i != N; ++i) {
        os << "state" << i;
        for (log_prob_dist::const_iterator iter = log_emission[i].begin();
             iter != log_emission[i].end(); ++iter) {
            os << "\t" << ((*iter >= cutoff) ? 0.0 : exp(*iter));
        }
        os << endl;
    }
    
    os.precision(logP_pre);
    os << "\nOverall log likelihood:\n" << logP << endl;
    
    os.flags(original_flgs);
    os.precision(original_pre);
}

NonparametricRatioHMM::log_likelihood
NonparametricRatioHMM::trainParas(const std::vector<obs_seq> &obs, const prob_dist &init_pi,
                                  const std::vector<prob_dist> &init_trans,
                                  const std::vector<prob_dist> &init_emission,
                                  bool verbose, std::ostream &os) {
    state = 0;
    size_t n_obs = 0;
    for (vector<obs_seq>::const_iterator iter = obs.begin();
         iter != obs.end(); ++iter) {
        if (iter -> size() == 0) {
            throw runtime_error("Any single observation sequence should contain at least 1 element!");
        }
        n_obs += iter -> size();
    }
    if (n_obs == obs.size()) {
        throw runtime_error("At least 1 observation sequence should contain at least 2 elements!");
    }
    setParas(init_pi, init_trans, init_emission);
    state = 0;
    
    const double cutoff = MARK + 2.0 * log(n_obs);
    switchMarks(MARK, 2.0 * cutoff);
    
    vector<vector<double> > log_binomial_coeff;
    log_choose(obs, log_binomial_coeff);
    
    vector<double> log_ratio(nbins + 1);
    log_ratio[0] = 2.0 * cutoff;
    for (size_t i = 1; i <= nbins; ++i) {
        log_ratio[i] = log(static_cast<double>(i) / nbins);
    }
    
    // allocate spaces for storing intermediate variables during the iterations.
    vector<Matrix> site_log_emission(obs.size()), alpha(obs.size()), gama(obs.size());
    for (size_t i = 0; i != obs.size(); ++i) {
        site_log_emission[i].resize(obs[i].size(), vector<double>(N, 0.0));
        alpha[i].resize(obs[i].size(), vector<double>(N, 0.0));
        gama[i].resize(obs[i].size(), vector<double>(N, 0.0));
    }
    log_prob_dist n_log_pi(log_pi);
    vector<log_prob_dist> n_log_trans(log_trans);
    vector<log_prob_dist> n_log_emission(log_emission);
    
    // iterations for training parameters.
    get_site_log_emission(site_log_emission, obs, log_binomial_coeff, log_ratio, cutoff);
    forward(alpha, site_log_emission, cutoff);
    log_likelihood logProb = getLogLikelihoodUtil(alpha, cutoff);
    if (logProb >= cutoff) {
        throw runtime_error("Given the parameters for the initialization of HMM, \
it's impossible to generate the observation sequences!");
    }
    
    bool converge = (logProb >= 0.0) ? true : false;
    for (size_t iter_n = 0; iter_n <= max_iterations; ++iter_n) {
        if (verbose) printInfo(iter_n, logProb, os, cutoff);
        if (converge) {
            if (verbose) os << "\nConverged!" << endl;
            break;
        }
        if (iter_n == max_iterations) {
            if (verbose) os << "\nReach maximum number of iterations!" << endl;
            break;
        }
        backward(gama, n_log_trans, n_log_pi, alpha, site_log_emission, cutoff);
        adjustParas(n_log_pi, cutoff, 2.0 * cutoff);
        for (size_t i = 0; i != N; ++i) adjustParas(n_log_trans[i], cutoff, 2.0 * cutoff);
        fit_emission(n_log_emission, gama, site_log_emission, obs, log_binomial_coeff, log_ratio, cutoff);
        for (size_t i = 0; i != N; ++i) adjustParas(n_log_emission[i], cutoff, 2.0 * cutoff);
        log_pi = n_log_pi; log_trans = n_log_trans; log_emission = n_log_emission;
        
        // end of a single re-estimation procedure.
        // calculate the overall likelihood based on newly trained parameters.
        get_site_log_emission(site_log_emission, obs, log_binomial_coeff, log_ratio, cutoff);
        forward(alpha, site_log_emission, cutoff);
        log_likelihood temp = getLogLikelihoodUtil(alpha, cutoff);
        
        // Sometimes, the theoretical increasing is violated due to a limited computational precision.
        // if (temp >= cutoff || temp < logProb) {
        //     throw runtime_error("Internal bugs occur! The overall likelihood isn't increasing!");
        // }
        if (temp >= cutoff) throw runtime_error("Internal bugs occur!");
        
        if (temp >= 0.0 ||
            std::abs(temp - logProb) < std::abs(logProb) * log_likelihood_tol) converge = true;
        logProb = temp;
    }
    switchMarks(cutoff, 2.0 * MARK);
    state = 1;
    return logProb;
}

NonparametricRatioHMM::log_likelihood
NonparametricRatioHMM::getLogLikelihood(const std::vector<obs_seq> &obs) const {
    if (state == 0) {
        throw runtime_error("HMM parameters haven't been set or trained yet!");
    }
    
    for (vector<obs_seq>::const_iterator iter = obs.begin();
         iter != obs.end(); ++iter) {
        if (iter -> size() == 0) {
            throw runtime_error("Any single observation sequence should contain at least 1 element!");
        }
    }
    
    const double cutoff = MARK;
    vector<vector<double> > log_binomial_coeff;
    log_choose(obs, log_binomial_coeff);
    
    vector<double> log_ratio(nbins + 1);
    log_ratio[0] = 2.0 * cutoff;
    for (size_t i = 1; i <= nbins; ++i) {
        log_ratio[i] = log(static_cast<double>(i) / nbins);
    }
    
    vector<Matrix> site_log_emission(obs.size()), alpha(obs.size());
    for (size_t i = 0; i != obs.size(); ++i) {
        site_log_emission[i].resize(obs[i].size(), vector<double>(N, 0.0));
        alpha[i].resize(obs[i].size(), vector<double>(N, 0.0));
    }
    
    get_site_log_emission(site_log_emission, obs, log_binomial_coeff, log_ratio, cutoff);
    forward(alpha, site_log_emission, cutoff);
    return getLogLikelihoodUtil(alpha, cutoff);
}

inline bool NonparametricRatioHMM::is_log0(log_likelihood x) {
    return x >= MARK;
}

NonparametricRatioHMM::log_likelihood
NonparametricRatioHMM::posterior_prob(const std::vector<obs_seq> &obs,
                                      std::vector<std::vector<prob_dist> > &post_dist_states) const {
    if (state == 0) {
        throw runtime_error("HMM parameters haven't been set or trained yet!");
    }
    
    size_t n_obs = 0;
    for (vector<obs_seq>::const_iterator iter = obs.begin();
         iter != obs.end(); ++iter) {
        if (iter -> size() == 0) {
            throw runtime_error("Any single observation sequence should contain at least 1 element!");
        }
        n_obs += iter -> size();
    }
    if (n_obs == 0) {
        throw runtime_error("No observations at all! Can't perform the posterior decoding procedure.");
    }
    
    const double cutoff = MARK;
    
    vector<vector<double> > log_binomial_coeff;
    log_choose(obs, log_binomial_coeff);
    
    vector<double> log_ratio(nbins + 1);
    log_ratio[0] = 2.0 * cutoff;
    for (size_t i = 1; i <= nbins; ++i) {
        log_ratio[i] = log(static_cast<double>(i) / nbins);
    }
    
    vector<Matrix> site_log_emission(obs.size()), alpha(obs.size()), gama(obs.size());
    post_dist_states.resize(obs.size());
    for (size_t i = 0; i != obs.size(); ++i) {
        site_log_emission[i].resize(obs[i].size(), vector<double>(N, 0.0));
        alpha[i].resize(obs[i].size(), vector<double>(N, 0.0));
        gama[i].resize(obs[i].size(), vector<double>(N, 0.0));
        post_dist_states[i].resize(obs[i].size(), vector<double>(N, 0.0));
    }
    log_prob_dist n_log_pi(log_pi);
    vector<log_prob_dist> n_log_trans(log_trans);
    
    // forward-backward procedures.
    get_site_log_emission(site_log_emission, obs, log_binomial_coeff, log_ratio, cutoff);
    forward(alpha, site_log_emission, cutoff);
    log_likelihood logProb = getLogLikelihoodUtil(alpha, cutoff);
    if (logProb >= cutoff) {
        throw runtime_error("Based on the HMM parameters, the overall likelihood of the \
observation sequences is 0! Can't perform the posterior decoding procedure.");
    }
    backward(gama, n_log_trans, n_log_pi, alpha, site_log_emission, cutoff);
    
    // filling in the posterior probability distributions.
    for (size_t i = 0; i != obs.size(); ++i) {
        Matrix::iterator post_iter = post_dist_states[i].begin(),
                         gama_iter = gama[i].begin();
        while (gama_iter != gama[i].end()) {
            for (size_t j = 0; j != N; ++j) {
                if (gama_iter -> at(j) < cutoff) post_iter -> at(j) = exp(gama_iter -> at(j));
            }
            ++post_iter;
            ++gama_iter;
        }
    }
    
    return logProb;
}

