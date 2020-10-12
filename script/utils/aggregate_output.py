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
import random
import numpy as np
import json

# project files
import utils

def callback(data, ecc_code, raw_data):
    if data['uid_ignoring_fields'] not in raw_data:
        raw_data[data['uid_ignoring_fields']] = []
    raw_data[data['uid_ignoring_fields']].append(data)

def aggregate_einsim_results(infiles, outfile, fields_to_combine, print_stdout, graph, bootstrap):
    if not outfile and not print_stdout and not graph:
        print("[ERROR] must either output to file, stdout, or graph!!")
        sys.exit(-1)

    # check if the user specified an output file. if not, stdout will be used
    if outfile and os.path.isfile(outfile):
        print("[ERROR] outfile already exists! will not crush it")
        sys.exit(-1)
    
    # parse the input files
    raw_data = {}
    print("Parsing", len(infiles), "files...")
    # raw_data = utils.parse_all_files(infiles, experimental=False)
    ecc_codes = utils.parse_files_incremental(infiles, False, fields_to_combine, callback, raw_data)
    ecc_codes_by_uid = {ecc_codes[code]['uid'] : code for code in ecc_codes}

    # coalesce the input data
    f = open(outfile, 'w') if outfile else None
    for uid_if, data in raw_data.items():
        combined_data = utils.combine_runs(data, fields_to_combine)
        out_str = utils.get_output_string(combined_data)
        if f:            f.write(out_str)
        if print_stdout: sys.stdout.write(out_str)
    if f: f.close()

    # can plot the output data if desired
    if graph:
        import plot_einsim_results
        plot_einsim_results.plot(raw_data)


    # bootstrap the data and generate CIs
    if bootstrap:
        # print the ECC code
        assert len(ecc_codes) == 1, "can only bootstrap for one ECC code!"
        print("[ECC]", json.dumps(list(ecc_codes.values())[0]))

        # print the statistic
        combined_statistic = combined_data['observations']['PER_BIT_ERROR_COUNT']['hist_databurst']
        print("[DATA] COMBINED:", np.array(combined_statistic) / sum(combined_statistic))
        
        # compute the CIs
        for uid_if, data in raw_data.items():
            nbs = int(1e9)
            print("[INFO] bootstrapping", nbs, "rounds of nsamples:", bootstrap, "from uid:", uid_if, "with", len(data), "data points to generate CIs") 
            assert bootstrap <= len(data), "cannot take more bootstrap samples than there are data points"
            ndata = len(data)
            for i in range(bootstrap): 
                # take a bootstrap sample - subsample ``bootstrap'' values from the dataset
                databurst_errors = np.array([0 for _ in data[0]['observations']['PER_BIT_ERROR_COUNT']['hist_databurst']], dtype=int)
                codeburst_errors = np.array([0 for _ in data[0]['observations']['PER_BIT_ERROR_COUNT']['hist_codeburst']], dtype=int)

                # determine how many samples fall in the nonzero (i.e., data) range - rest are 0, so they don't need any work
                n_in_range = np.random.binomial(nbs, ndata / float(nbs))
                print("[DEBUG] drawing", n_in_range, "samples from the data: {0:.2f}".format(100.0 * i / float(bootstrap - 1)) + '% complete')
                for s in range(n_in_range):
                    sample_idx = random.randint(0, len(data) - 1)
                    databurst_errors += data[sample_idx]['observations']['PER_BIT_ERROR_COUNT']['hist_databurst']
                    codeburst_errors += data[sample_idx]['observations']['PER_BIT_ERROR_COUNT']['hist_codeburst']
 
                # compute the statistic
                print('[BOOTSTRAP]', (list(databurst_errors), list(codeburst_errors)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply the EIN inference methodology using the results of EINSim')
    parser.add_argument(dest='input_filenames', type=str, nargs='*',
                        help='input data files to parse (i.e., output file(s) of EINSim)')
    parser.add_argument('-o', '--outfile', type=str, default=None,
                        help='output filename for aggregated input data')
    parser.add_argument('-d', '--dirname', type=str, default=None,
                        help='local directory containing input files (scanned recursively)')
    parser.add_argument('-b', '--bootstrap', type=int, default=None,
                        help='number of boostrap samples to create confidence intervals (requires single-sample outputs)')
    parser.add_argument('-s', '--suppress-output', action='store_true', default=False,
                        help='suppress output to stdout (will still go to file if requested)')
    parser.add_argument('-c', '--combine', type=str, action='append',
                        help='fields to combine when analyzing multiple dumps')
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
    if not args.combine:
        args.combine = []
    aggregate_einsim_results(clean_filenames, args.outfile, args.combine, not args.suppress_output, args.graph, args.bootstrap)
