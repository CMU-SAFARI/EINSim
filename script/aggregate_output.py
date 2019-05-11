#!/usr/bin/env python3
# -*- coding:utf-8 -*-
###
# @file aggregate_output.py
#
# @brief Aggregates multiple EINSim outputs into one coalesced file
#
# @author Minesh Patel
# Contact: minesh.patelh@gmail.com
import sys
import os
import argparse

# project files
import utils

def graph_output(raw_data):
    import warnings
    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    print("Displaying plot of data")
    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    for run_uid, components in raw_data.items():
        ecc_category, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp = components['metadata']

        histogram = components['re_data']
        x, y = zip(*sorted(histogram.items()))
        y = np.array(y) / np.sum(y)
        ax.plot(x, y, label='re_' + run_uid, color='r')

        histogram = components['ue_data']
        x, y = zip(*sorted(histogram.items()))
        y = np.array(y) / np.sum(y)
        ax.plot(x, y, label='ue_' + run_uid, color='g')
    ax.set_xlabel('# Errors per Word')
    ax.set_ylabel('Relative Frequency')
    ax.set_yscale('log')
    fig.tight_layout()
    plt.show()

def condense_raw_data(raw_data):
    print("Condensing raw data...")

    all_ecc_types = {}
    for run_uid, components in raw_data.items():
        ecc_category, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp = components['metadata']
        uber = components['uber']
        rber = components['rber']

        ecc_id = (ecc_category, t, k, n, m, bl, bcl, ps, ed, cd, dp)
        if ecc_id not in all_ecc_types:
            all_ecc_types[ecc_id] = []
        all_ecc_types[ecc_id].append((rber, uber))

    return all_ecc_types

def aggregate_einsim_results(infiles, outfile, print_stdout, graph):
    if not outfile and not print_stdout and not graph:
        print("[ERROR] must either output to file, stdout, or graph!!")
        sys.exit(-1)

    # check if the user specified an output file. if not, stdout will be used
    if outfile and os.path.isfile(outfile):
        print("[ERROR] outfile already exists! will not crush it")
        sys.exit(-1)
    
    # parse the input files
    print("Parsing", len(infiles), "files...")
    raw_data = utils.parse_all_files(infiles, experimental=False)

    # coalesce the input data
    out_str = utils.get_minimized_output_string(raw_data)
    if outfile:
        with open(outfile, 'w') as f:
            f.write(out_str)
    if print_stdout:
        sys.stdout.write(out_str)

    # can plot the output data if desired
    if graph:
        graph_output(raw_data)

    # Debug condensed output mode
    if False:
        condensed_data = condense_raw_data(raw_data)
        for ecc_type in condensed_data:
            rbers, ubers = zip(*condensed_data[ecc_type])
            print(ecc_type, rbers, ubers)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply the EIN inference methodology using the results of EINSim')
    parser.add_argument(dest='input_filenames', type=str, nargs='*',
                        help='input data files to parse (i.e., output file(s) of EINSim)')
    parser.add_argument('-o', '--outfile', type=str, default=None,
                        help='output filename for aggregated input data')
    parser.add_argument('-d', '--dirname', type=str, default=None,
                        help='local directory containing input files (scanned recursively)')
    parser.add_argument('-s', '--suppress-output', action='store_true', default=False,
                        help='suppress output to stdout (will still go to file if requested)')
    parser.add_argument('-g', '--graph', action='store_true', default=False,
                        help='enable graphing')
    args = parser.parse_args(sys.argv[1:])
    
    # list of all input files
    all_filenames = args.input_filenames

    # if a directory was specified, recurse into it and add all files to the list
    if args.dirname:
        for dp, dn, fnames in os.walk(args.dirname):
            for f in fnames:
                all_filenames.append(os.path.join(dp, f))
    if len(all_filenames) == 0:
        print("[ERROR] must have at least one input file to parse")
        sys.exit(-1)
    print("[INFO] parsing", len(all_filenames), "input files")

    # test all input files for existence (sanity)
    clean_filenames = []
    for fname in all_filenames:
        if not os.path.isfile(fname):
            print("[ERROR] invalid input filename: \"" + fname + "\"")
            sys.exit(-1)
        else:
            clean_filenames.append(fname)
 
    # reduce and output the results as required
    aggregate_einsim_results(clean_filenames, args.outfile, not args.suppress_output, args.graph)
