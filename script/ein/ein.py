#!/usr/bin/env python3
# -*- coding:utf-8 -*-
###
# @file ein.py
#
# @brief Applies EINSim output to real experimental data
#
# @author Minesh Patel
# Contact: minesh.patelh@gmail.com
import sys
import os
import argparse
import subprocess
import pipes
import getpass
import time
import traceback
import math
import heapq
import bisect
import random
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import numpy as np
from scipy.stats import multinomial

# project files
import utils

log_probabilities = {}
max_n_models_to_plot = 20
top_n_models_heapq = {}
top_n_models_data = {}

def callback(sim_uid, sim_data, *args):
    #
    # the experiments 
    #
    # for arg in args:
    #     print ("LEN:", arg)
    exp_uid, exp_data, o_sample_histograms = args
    n_bootstrap_samples = len(o_sample_histograms)
    o_ecc_category, o_t, o_k, o_n, o_m, o_rber, o_bl, o_bcl, o_ps, o_ed, o_cd, o_dp = exp_data['metadata']

    #
    # the simulation 
    #
    s_ecc_category, s_t, s_k, s_n, s_m, s_rber, s_bl, s_bcl, s_ps, s_ed, s_cd, s_dp = sim_data['metadata']
    s_uber = sim_data['uber']
    s_rber = sim_data['rber']
    if o_bl != s_bl:
        print("[WARNING] burst length mismatch! Skipping model", exp_uid, sim_uid)
        return

    s_xlabel = s_ecc_category + '(' + str(s_n) + ', ' + str(s_k) + ', ' + str(s_t) + ') [r:' + ('%s' % float('%.4g' % s_rber)) + ' u:' + ('%s' % float('%.4g' % s_uber)) + ']'
    s_histogram = sim_data['ue_data']

    # obtain the probabilities of the simulation - these are the single-trial probabilities
    s_histogram_sum = float(sum(s_histogram.values()))
    single_trial_probabilities = {n_errs : count / s_histogram_sum for n_errs, count in s_histogram.items()}
    for i in range(s_bl + 1):
        if i not in single_trial_probabilities:
            single_trial_probabilities[i] = 1e-10 # for numerical stability when computing the multinomial coefficients

    # compute a multinomial PMF
    s_pmf = [prob for n_errs, prob in sorted(single_trial_probabilities.items())]

    # run N bootstrap samples as required
    log_probs = []
    n_iterations = 1 if n_bootstrap_samples == 0 else n_bootstrap_samples # bootstrapping does NOT use the original sample
    for o_sample_histogram in o_sample_histograms: 
        # obtain the probabilities of the observations - this is for input to a multinomial distribution
        o_sample_histogram_sum = sum(o_sample_histogram.values())
        # o_sample_probs = {n_errs : count / o_sample_histogram_sum for n_errs, count in o_sample_histogram.items()}
        o_sample_probs = {n_errs : count for n_errs, count in o_sample_histogram.items()}
        for i in range(o_bl + 1):
            if i not in o_sample_probs:
                o_sample_probs[i] = 1e-10 # for numerical stability when computing the multinomial coefficients

        o_sample_pmf = [prob  for n_errs, prob in sorted(o_sample_probs.items())]

        x = np.array(o_sample_pmf, dtype=int)
        n = int(o_sample_histogram_sum)
        p = np.array(s_pmf, dtype=float)            

        # WARNING
        # there is a BUG IN SCIPY MULTINOMIAL that results in NAN being generated
        #   when any of the PMF probabilities are 0
        # therefore, we skip all the nan/inf validation steps and compute the PMF directly
        from scipy.special import xlogy, gammaln
        log_prob = gammaln(n+1) + np.sum(xlogy(x, p) - gammaln(x+1), axis=-1)
        # log_prob = multinomial.logpmf(o_sample_pmf, n=int(o_sample_histogram_sum), p=s_pmf)
        # print prob, log_prob
        # log_prob = float(log_prob) if not math.isnan(log_prob) else float('-inf')
        log_probs.append(log_prob)
    
    # keep the log-probability details for all models
    global log_probabilities
    if exp_uid not in log_probabilities:
        log_probabilities[exp_uid] = {}
    assert sim_uid not in log_probabilities[exp_uid], "Repeat simulation: " + sim_uid + " == " + str(sim_data)
    log_probabilities[exp_uid][sim_uid] = [log_probs, s_xlabel]
    
    # keep the plot details for the top 20 models
    global top_n_models_data
    global top_n_models_heapq
    if exp_uid not in top_n_models_data:
        top_n_models_data[exp_uid] = {"exp_data" : exp_data, "sims" : {}}
        top_n_models_heapq[exp_uid] = []
    sim_max_logprob = np.max(log_probs)
    if len(top_n_models_heapq[exp_uid]) < max_n_models_to_plot or top_n_models_heapq[exp_uid][0][0] < sim_max_logprob:
        if len(top_n_models_heapq[exp_uid]) == max_n_models_to_plot:
            heapitem_logprob, heapitem_uid = heapq.heappop(top_n_models_heapq[exp_uid])
            top_n_models_data[exp_uid]["sims"].pop(heapitem_uid)
        heapq.heappush(top_n_models_heapq[exp_uid], (sim_max_logprob, sim_uid))
        top_n_models_data[exp_uid]["sims"][sim_uid] = (log_probs, s_xlabel, sim_data)


def evaluate_likelihoods_incremental(observations, simulation_filenames, n_bootstrap_samples, graph):
    if n_bootstrap_samples == 0:
        print("[INFO] Bootstrapping DISABLED")
    else:
        print("[INFO] Bootstrapping with", n_bootstrap_samples, "samples")
        assert n_bootstrap_samples >= 1, "[ERROR] number of bootstrap samples must be positive"

    for observation_uid, observation_components in observations.items():
        o_ecc_category, o_t, o_k, o_n, o_m, o_rber, o_bl, o_bcl, o_ps, o_ed, o_cd, o_dp = observation_components['metadata']
        o_uber = observation_components['uber']
        o_histogram = observation_components['ue_data']

        # resample the histogram to create all the bootstrap variants
        o_sample_histograms = []
        if n_bootstrap_samples == 0:
            o_sample_histograms.append(o_histogram)
        else:
            # create a PMF of the observations
            hist_sum = np.sum(o_histogram.values())
            categories, pmf = zip(*[(n_errs, count) for n_errs, count in sorted(o_histogram.items()) if count != 0])
            pmf_sum = float(np.sum(pmf))
            pmf_normed = np.array(pmf) / pmf_sum
            rv = multinomial(n=pmf_sum, p=pmf_normed)
            samples = rv.rvs(n_bootstrap_samples)
            # print pmf_normed
            # print sample
            for bs_sample in samples: 
                new_hist = {n_errs : 0 for n_errs in o_histogram.keys()}
                for n_errs, count in zip(categories, bs_sample):
                    new_hist[n_errs] = count
                o_sample_histograms.append(new_hist)
        # print o_sample_histograms
        # sys.exit(0)

        # compute for all probabilities against all runs and NORMALIZE (note: BL must the same!)
        print("Computing probabilities with", n_bootstrap_samples, "bootstrap samples for observations:", observation_uid)
        utils.parse_files_incremental(simulation_filenames, False, callback, observation_uid, observation_components, o_sample_histograms)

# use log-probabilities to evaluate this 
def evaluate_likelihoods(observations, simulations, n_bootstrap_samples, graph):
    if n_bootstrap_samples == 0:
        print("[INFO] Bootstrapping DISABLED")
    else:
        print("[INFO] Bootstrapping with", n_bootstrap_samples, "samples")
        assert n_bootstrap_samples >= 1, "[ERROR] number of bootstrap samples must be positive"

    for observation_uid, observation_components in observations.items():
        o_ecc_category, o_t, o_k, o_n, o_m, o_rber, o_bl, o_bcl, o_ps, o_ed, o_cd, o_dp = observation_components['metadata']
        o_uber = observation_components['uber']
        o_histogram = observation_components['ue_data']

        # resample the histogram to create all the bootstrap variants
        o_sample_histograms = []
        if n_bootstrap_samples == 0:
            o_sample_histograms.append(o_histogram)
        else:
            # create a PMF of the observations
            hist_sum = np.sum(o_histogram.values())
            categories, pmf = zip(*[(n_errs, count) for n_errs, count in sorted(o_histogram.items()) if count != 0])
            pmf_sum = float(np.sum(pmf))
            pmf_normed = np.array(pmf) / pmf_sum
            rv = multinomial(n=pmf_sum, p=pmf_normed)
            samples = rv.rvs(n_bootstrap_samples)
            # print pmf_normed
            # print sample
            for bs_sample in samples: 
                new_hist = {n_errs : 0 for n_errs in o_histogram.keys()}
                for n_errs, count in zip(categories, bs_sample):
                    new_hist[n_errs] = count
                o_sample_histograms.append(new_hist)
        # print o_sample_histograms
        # sys.exit(0)

        # compute for all probabilities against all runs and NORMALIZE (note: BL must the same!)
        print("Computing probabilities over", len(simulations), "simulations with", n_bootstrap_samples, "bootstrap samples for observations:", observation_uid)
        for simulation_run_uid, simulation_components in simulations.items():
            callback(simulation_run_uid, simulation_components, observation_uid, observation_components, o_sample_histograms)

# plot the N-best distributions
# we generate TWO plots - the log-probabilities and the distributions
def plot_likelihoods(exp_uid, exp_data, sims):    
    # expand the experimental data
    e_ecc_category, e_t, e_k, e_n, e_m, e_rber, e_bl, e_bcl, e_ps, e_ed, e_cd, e_dp = exp_data['metadata']
    e_uber = exp_data['uber']
    e_uber_hist = exp_data['ue_data']
    e_rber = exp_data['rber']
    e_rber_hist = exp_data['re_data']
    
    #
    # show a bar-chart of the models' log-probabilities
    #
    ratio = 4
    fig, ax = plt.subplots(1, 1, figsize=(5, 7))
    sim_uids, sim_metadatas = zip(*sorted(sims.items(), key=lambda x: np.median(x[1][0]))) 
    log_probs, s_xlabels, sim_datas = zip(*sim_metadatas)
    print(log_probs)

    yvals = list(np.negative(log_probs))
    ymeds = np.median(yvals, axis=1)
    ymaxs = np.max(yvals, axis=1)
    ymins = np.min(yvals, axis=1)
    xvals = np.arange(len(yvals))
    bar_width = 0.35
    ax.bar(xvals, ymeds, yerr=(ymeds - ymins, ymaxs - ymeds), width=bar_width, ecolor='black', capsize=4)
    ax.set_xticklabels(s_xlabels)
    ax.set_ylabel('-Log(Likelihood)')
    def exp_notation_formatter(x, p):
        return ('%.1f' % (x / (10 ** 7))) + 'e7'
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(exp_notation_formatter))
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(exp_notation_formatter))
    ax.tick_params(axis='x', labelrotation=90)
    fig.text(0.6, 0.03, 'ECC Schemes', ha='center')
    fig.tight_layout()
    fig.canvas.set_window_title(exp_uid)

    #
    # show the lineplot of the models' actual histogram data
    #
    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    x = range(e_bl + 1)
    for s_uid, log_probs, s_xlabel, sim_data in zip(*[sim_uids, log_probs, s_xlabels, sim_datas]):
        if s_uid == sim_uids[-1]: # assuming -1 is the best-fit
            lw = 2
            best_fit = True
            ulabel = 'simulated post-' + s_xlabel
            rlabel = 'simulated pre-' + s_xlabel
            col = plt.rcParams['axes.prop_cycle'].by_key()['color'][1]

            # plot second-to-last so z-ordering is correct
            x, y = zip(*sorted(e_uber_hist.items()))
            label = 'Experimental Data'
            # label = 'observed ' + ('[%s]' % float('%.3g' % o_uber))
            ax.plot(x, np.array(y) / float(sum(y)), ls='-', lw=2, label=label)

            ulabel = 'Maximum-Likelihood Model'
        else:
            lw = 0.5
            best_fit = False
            ulabel = None
            rlabel = None
            col = 'grey'

        # plot UBER
        x, y = zip(*sorted(sim_data["ue_data"].items()))
        p = ax.plot(x, np.array(y) / float(sum(y)), ls=':', lw=lw, label=ulabel, color=col)
        color = p[0].get_color()

        # plot RBER
        # if best_fit:
        #     x, y = zip(*sorted(sim_data["re_data"].items()))
        #     ax.plot(x, np.array(y) / float(sum(y)), color=color, ls='--', lw=lw, label=rlabel)

    ax.set_xlabel('Number of Bit Errors in a ' + str(e_bl) + '-bit Word')
    ax.set_ylabel('Probability')
    ax.set_yscale('log')
    ax.plot([], [], ls=':', lw=0.5, label='Lower-Likelihood Models', color='grey')
    ax.legend()
    fig.tight_layout()
    fig.canvas.set_window_title(exp_uid)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply the EIN inference methodology using the results of EINSim')
    parser.add_argument(dest='simulation_datafile', type=str, nargs='*',
                        help='input data files to parse (i.e., output file(s) of EINSim)')
    parser.add_argument('-d', '--simulation_dirname', type=str, default=None,
                        help='directory containing simulation files (scanned recursively)')
    parser.add_argument('-o', '--outfile', type=str, default=None,
                        help='output filename for evaluated observations')
    parser.add_argument('-g', '--graph', action='store_true', default=False,
                        help='enable graphing')
    parser.add_argument('-e', '--experimental_datafile', type=str, action='append', default=[],
                        help='filename containing the experimental data to evaluate (can provide more than once)')
    parser.add_argument('-b', '--bootstrap', type=int, default=0,
                        help='number of bootstrap samples to evaluate confidence intervals with (local use only)')
    parser.add_argument('-c', '--condense', action='store_true', default=False,
                        help='determine the RBER -> UBER relationship in a condensed form (local use only)')
    parser.add_argument('-s', '--suppress-output', action='store_true', default=False,
                        help='suppress output to stdout (will still go to file if requested)')
    args = parser.parse_args(sys.argv[1:])

    # check if the user specified an output file. if not, stdout will be used
    if args.outfile and isfile(args.outfile):
        print("[ERROR] outfile already exists! will not crush it")
        sys.exit(-1)

    # gather the simulation filenames
    simulation_filenames = args.simulation_datafile
    if args.simulation_dirname:
        for dp, dn, fnames in os.walk(args.simulation_dirname):
            for f in fnames:
                simulation_filenames.append(os.path.join(dp, f))
    
    if len(simulation_filenames) == 0:
        print("[ERROR] must have at least one input file to parse")
        parser.print_help()
        sys.exit(-1)

    clean_sim_filenames = []
    for fname in simulation_filenames:
        if not os.path.isfile(fname):
            print("[ERROR] invalid input filename: \"" + fname + "\"")
            sys.exit(-1)
        else:
            clean_sim_filenames.append(fname)

    # gather the experiment filenames
    experimental_filenames = args.experimental_datafile
    
    if len(experimental_filenames) == 0:
        print("[ERROR] must have at least one experimental file to parse")
        parser.print_help()
        sys.exit(-1)
    
    clean_exp_filenames = []
    for fname in experimental_filenames:
        if not os.path.isfile(fname):
            print("[ERROR] invalid input filename: \"" + fname + "\"")
            sys.exit(-1)
        else:
            clean_exp_filenames.append(fname)

    # parse the experimental input files
    print("Parsing experimental results", len(clean_exp_filenames), "files...")
    raw_data_exp = utils.parse_all_files(clean_exp_filenames, experimental=True)
 
    # we default to the 'incremental' pass mode because it works for large input datasets
    # however, the normal pass that reads the entire input sets into memory also works
    if True:
        evaluate_likelihoods_incremental(raw_data_exp, clean_sim_filenames, args.bootstrap, args.graph)
    else:
        # parse the simulation input files
        print("Parsing simulation results", len(clean_sim_filenames), "files...")
        raw_data_sim = utils.parse_all_files(clean_sim_filenames, experimental=False)

        # evaluate EIN over the input experimental data
        evaluate_likelihoods(raw_data_exp, raw_data_sim, args.bootstrap, args.graph)

    # print the probabilities if requested
    if not args.suppress_output:
        for exp_uid in log_probabilities:
            print("[INFO] Evaluations for experiment:", exp_uid)
            sim_uids, _ = zip(*sorted(log_probabilities[exp_uid].items(), key=lambda x: np.median(x[1][0]))[-100:])
            for sim_uid in sim_uids:
                log_probs, s_xlabel = log_probabilities[exp_uid][sim_uid]
                print("[DATA] log(p): [", min(log_probs), max(log_probs), "] " + s_xlabel + " -> UID: " + sim_uid)
            
    # graph the probabilities if requested
    if args.graph:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib import ticker
        from matplotlib.patches import Patch

        for exp_uid in top_n_models_data:
            exp_data = top_n_models_data[exp_uid]["exp_data"]
            sims = top_n_models_data[exp_uid]["sims"]
            plot_likelihoods(exp_uid, exp_data, sims)
        plt.show()
    sys.exit(0)
