#!/usr/bin/env python3
# -*- coding:utf-8 -*-
###
# @file utils.py
#
# @brief Common helper functions
#
# @author Minesh Patel
# Contact: minesh.patelh@gmail.com
import sys
import os
import copy
import traceback
import json
try:
    import lzma
except ImportError:
    print("[WARNING] unable to import lzma module required to read lzma/xz files. Use Python3.3+ or backports lzma for this functionality")
import gzip
import bz2


OBS_N_ERRORS_PER_BURST = 'N_ERRORS_PER_BURST'
OBS_PER_BIT_ERROR_COUNT = 'PER_BIT_ERROR_COUNT'
OBSERVABLES = {OBS_N_ERRORS_PER_BURST, OBS_PER_BIT_ERROR_COUNT}

# enum to index the metadata
MI_UID = 0
MI_NW = 1
MI_BL = 2
MI_BCL = 3
MI_PS = 4
MI_EM = 5
MI_CD = 6
MI_DP = 7
MI_CDP = 8
MI_OBS = 9
MI_MAX = 10

def add_to_hist(base, addend):
    for before, after, count in addend:
        if before not in base:
            base[before] = {}
        if after not in base[before]:
            base[before][after] = 0
        base[before][after] += count

def open_sim_file(fname):
    def warn_about_compressed_progress_bar():
        print("\t[WARNING] progress bar will go >100% for compressed file. On the bright side, the higher it goes, the better your compression ratio is :D")
    
    if fname[-3:] == '.xz':
        warn_about_compressed_progress_bar()
        return lzma.open(fname, mode='rt')
    elif fname[-3:] == '.gz':
        warn_about_compressed_progress_bar()
        return gzip.open(fname, mode='rt')
    elif fname[-4:] == '.bz2':
        warn_about_compressed_progress_bar()
        return bz2.open(fname, mode='rt')
    else:
        return open(fname, 'r')

# example input lines
# [DATA] REP: p:584576191 t:1 k:4 n:12 rber:1e-06 bl:256 bcl:768 ps:0 ed:UNIFORM_RANDOM cd:ALL_TRUE_OR_ALL_ANTI dp:RANDOM [ 0:9998:10000 1:2:0 ]
# [DATA] REP: p:584576191 t:3 k:4 n:28 rber:0.0331251 bl:256 bcl:1792 ps:0 ed:UNIFORM_RANDOM cd:ALL_TRUE_OR_ALL_ANTI dp:RANDOM [ 0:0:9051 1:0:900 2:0:49 31:2:0 32:2:0 33:1:0 34:3:0 35:4:0 36:9:0 37:16:0 38:16:0 39:24:0 40:34:0 41:44:0 42:70:0 43:72:0 44:98:0 45:131:0 46:135:0 47:182:0 48:190:0 49:234:0 50:252:0 51:309:0 52:276:0 53:365:0 54:407:0 55:413:0 56:426:0 57:437:0 58:449:0 59:431:0 60:499:0 61:458:0 62:426:0 63:357:0 64:402:0 65:353:0 66:365:0 67:303:0 68:298:0 69:252:0 70:203:0 71:173:0 72:150:0 73:146:0 74:114:0 75:98:0 76:79:0 77:59:0 78:55:0 79:44:0 80:25:0 81:31:0 82:16:0 83:18:0 84:15:0 85:7:0 86:4:0 87:5:0 88:2:0 89:3:0 92:4:0 93:1:0 95:1:0 99:1:0 101:1:0 ]
def parse_einsim_output_line(line, fields_to_combine, is_experimental=False):
    if len(line) > 6 and line[:6] == "[DATA]" in line:
        meta, raw_data = line[7:].strip().split('[')
        meta_split = meta.split(' ')
        metadata_extracted = dict([i.split(':', 1) for i in meta_split[0:-1]])

        metadata = [None for i in range(MI_MAX)]
        metadata[MI_UID] = -1 if 'uid' not in metadata_extracted else int(metadata_extracted['uid'])
        metadata[MI_NW ] = -1 if 'nw'  not in metadata_extracted else int(metadata_extracted['nw'])
        metadata[MI_BL ] = -1 if 'bl'  not in metadata_extracted else int(metadata_extracted['bl'])
        metadata[MI_BCL] = -1 if 'bcl' not in metadata_extracted else int(metadata_extracted['bcl'])
        metadata[MI_PS ] = -1 if 'ps'  not in metadata_extracted else int(metadata_extracted['ps'])
        metadata[MI_EM ] = -1 if 'em'  not in metadata_extracted else metadata_extracted['em']
        metadata[MI_CD ] = -1 if 'cd'  not in metadata_extracted else metadata_extracted['cd']
        metadata[MI_DP ] = -1 if 'dp'  not in metadata_extracted else metadata_extracted['dp']
        metadata[MI_CDP] = -1 if 'cdp' not in metadata_extracted else metadata_extracted['cdp']
        metadata[MI_OBS] = -1 if 'obs' not in metadata_extracted else metadata_extracted['obs']
        assert 'obs' in metadata_extracted, "No observational data for run"
        obs = metadata_extracted['obs']
        
        # ensure that experiments get unique UIDs - we do not aggregate them
        # metadata = [ecc_category, p, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp]
        run_uid = run_uid_from_metadata(metadata)
        run_uid_ignoring_fields = run_uid_from_metadata_ignoring_fields(metadata, fields_to_combine)
        cur_experiment_id = None
        if is_experimental:
            cur_experiment_id = experiment_id
            run_uid += '_experiment_' + str(cur_experiment_id)
            run_uid_ignoring_fields += '_experiment_' + str(cur_experiment_id)
            experiment_id += 1

        # print(ecc_type, permutation, nerrs_correctable, n_data_bits, n_code_bits, dp, error_model)
        # print(hist_data)
        
        if obs not in OBSERVABLES:
            print("[WARNING] unknown observable \"" + obs + "\"")
            return None 

        # assemble the final structure summarizing the run
        # omit the permutation - we're aggregating over it!
        parsed = {'metadata' : metadata
                   , 'experiment_id' : cur_experiment_id
                   , 'uid' : run_uid
                   , 'uid_ignoring_fields' : run_uid_ignoring_fields
                   , 'observations' : {}
        }


        # parse the different observable types and add them to the structure
        if obs == OBS_N_ERRORS_PER_BURST:
            data = [map(int, entry.split(':')) for entry in raw_data.split(' ') if ':' in entry]
            
            raw_histogram = {}
            re_histogram = {}
            ue_histogram = {}
            for n_errs, re_count, ue_count in data:
                if n_errs not in raw_histogram:
                    raw_histogram[n_errs] = [0, 0]
                raw_histogram[n_errs][0] += re_count
                raw_histogram[n_errs][1] += ue_count

                if re_count > 0:
                    assert n_errs <= bcl, "Cannot have more errors in a burst_codeword than there are bits: errs:" + str(n_errs) + ", bcl:" + str(bcl)
                    if n_errs not in re_histogram:
                        re_histogram[n_errs] = 0
                    re_histogram[n_errs] += re_count
                if ue_count > 0:
                    assert n_errs <= bl, "Cannot have more errors in a burst than there are bits: errs:" + str(n_errs) + ", bl:" + str(bl)
                    if n_errs not in ue_histogram:
                        ue_histogram[n_errs] = 0
                    ue_histogram[n_errs] += ue_count

            # compute the RBER and UBER
            if re_histogram:
                # account for the pad size in the RBER - a real device wouldn't actually use those bits!
                rber_calculated = sum([n_errs * count for n_errs, count in re_histogram.items()]) / float((bcl - ps) * sum(re_histogram.values()))
            else:
                rber_calculated = -1.0
            uber_calculated = sum([n_errs * count for n_errs, count in ue_histogram.items()]) / float(bl * sum(ue_histogram.values()))
        
            # add the observations
            parsed['observations'][OBS_N_ERRORS_PER_BURST] = {
                  'rber' : rber_calculated
                , 'uber' : uber_calculated
                , 're_data' : re_histogram
                , 'ue_data': ue_histogram
                , 'raw_data' : raw_histogram
            }

        elif obs == OBS_PER_BIT_ERROR_COUNT:
            raw_pre, raw_post = raw_data.split(':')
            bit_hist_databurst = [int(i) for i in raw_pre.split(' ') if i]
            bit_hist_codeburst = [int(i) for i in raw_post[:-1].split(' ') if i]
            
            # add the observations
            parsed['observations'][OBS_PER_BIT_ERROR_COUNT] = {
                  'hist_databurst' : bit_hist_databurst
                , 'hist_codeburst' : bit_hist_codeburst
            }

        return parsed
    else:
        return None

def combine_observation(observable, runs):
    if observable == OBS_N_ERRORS_PER_BURST:
        re_histogram = {}
        ue_histogram = {}
        raw_histogram = {}
        for run in runs:
            run_data = run['observations'][observable]
            for n_errs, re_count, ue_count in run_data['raw_data']:
                if n_errs not in raw_histogram:
                    raw_histogram[n_errs] = [0, 0]
                raw_histogram[n_errs][0] += re_count
                raw_histogram[n_errs][1] += ue_count

                if re_count > 0:
                    assert n_errs <= run['metadata'][MI_BCL], "Cannot have more errors in a burst_codeword than there are bits: errs:" + str(n_errs) + ", bcl:" + str(run['metadata'][MI_BCL])
                    if n_errs not in re_histogram:
                        re_histogram[n_errs] = 0
                    re_histogram[n_errs] += re_count
                if ue_count > 0:
                    assert n_errs <= run['metadata'][MI_BL], "Cannot have more errors in a burst than there are bits: errs:" + str(n_errs) + ", bl:" + str(run['metadata'][MI_BL])
                    if n_errs not in ue_histogram:
                        ue_histogram[n_errs] = 0
                    ue_histogram[n_errs] += ue_count
            
        # compute the RBER and UBER for the final combined histogram
        if re_histogram:
            # account for the pad size in the RBER - a real device wouldn't actually use those bits!
            rber_calculated = sum([n_errs * count for n_errs, count in re_histogram.items()]) / float((bcl - ps) * sum(re_histogram.values()))
        else:
            rber_calculated = -1.0
        uber_calculated = sum([n_errs * count for n_errs, count in ue_histogram.items()]) / float(bl * sum(ue_histogram.values()))
        

        return {
              'rber' : rber_calculated
            , 'uber' : uber_calculated
            , 're_data' : re_histogram
            , 'ue_data': ue_histogram
            , 'raw_data' : raw_histogram
        }
    elif observable == OBS_PER_BIT_ERROR_COUNT:
        bit_hist_databurst = []
        bit_hist_codeburst = []
        for run in runs:
            a = run['observations'][observable]
            if bit_hist_databurst or bit_hist_codeburst:
                assert len(bit_hist_databurst) == len(a['hist_databurst']), "Cannot combine different lengths: " + str(len(bit_hist_databurst)) + "vs" + str(len(a['hist_databurst']))
                assert len(bit_hist_codeburst) == len(a['hist_codeburst']), "Cannot combine different lengths: " + str(len(bit_hist_codeburst)) + "vs" + str(len(a['hist_codeburst']))
                for i in range(len(bit_hist_databurst)): bit_hist_databurst[i] += a['hist_databurst'][i]
                for i in range(len(bit_hist_codeburst)): bit_hist_codeburst[i] += a['hist_codeburst'][i]
            else:
                for i in range(len(a['hist_databurst'])): bit_hist_databurst.append(a['hist_databurst'][i])
                for i in range(len(a['hist_codeburst'])): bit_hist_codeburst.append(a['hist_codeburst'][i])
        return {
              'hist_databurst' : bit_hist_databurst
            , 'hist_codeburst' : bit_hist_codeburst
        }
    else:
        assert False, "Unimplemented observable:" + observable


def combine_runs(data, fields_to_combine):
    for run in data[1:]:
        assert data[0]['uid_ignoring_fields'] == run['uid_ignoring_fields'], str(data[0]['uid_ignoring_fields']) + " != " + str(run)

    combined_run = {
          'metadata' : [None for _ in range(MI_MAX)]
        , 'experiment_id' : -1
        , 'uid' : None
        , 'uid_ignoring_fields' : None
        , 'observations' : {}
    }

    def combine_metadata_element(data, elem, mixed_pass=None, mixed_warn=None, mixed_error=None):
        all_elem  = {run['metadata'][elem] for run in data}
        if len(all_elem) > 1: 
            if mixed_error != None:
                print("[ERROR] trying to combinine different types of " + elem + ":", all_elem)
                assert False, "[ERROR] trying to combinine different types of " + elem + ": " + str(all_elem)
            elif mixed_warn != None:
                print("[WARN] combining different types of " + elem + " - result will be marked as " + str(ifmixed))
                return mixed_warn
            elif mixed_pass != None:
                return mixed_pass
            else:
                assert False, "Must define either warn, error, or pass on mixed"
        else:
            return all_elem.pop() # okay to pop - we're not using this set anymore

    combined_run['metadata'][MI_UID] = combine_metadata_element(data, MI_UID, mixed_pass=-1)
    word_list = [run['metadata'][MI_NW] for run in data]    
    invalid_word_list = any([i == -1 for i in word_list])
    if not invalid_word_list:
        combined_run['metadata'][MI_NW]  = sum(word_list)
        total_bits_measured = 0
        for run in data:
            n_bits_measured = run['metadata'][MI_BCL] * run['metadata'][MI_NW]
            total_bits_measured += n_bits_measured
    else:
        combined_run['metadata'][MI_NW]    = combine_metadata_element(data, MI_NW, mixed_pass=-1)
    combined_run['metadata'][MI_BL]  = combine_metadata_element(data, MI_BL,  mixed_error=-1)
    combined_run['metadata'][MI_BCL] = combine_metadata_element(data, MI_BCL, mixed_warn=-1)
    combined_run['metadata'][MI_PS]  = combine_metadata_element(data, MI_PS,  mixed_warn=-1)
    combined_run['metadata'][MI_EM]  = combine_metadata_element(data, MI_EM,  mixed_pass=-1)
    combined_run['metadata'][MI_CD]  = combine_metadata_element(data, MI_CD,  mixed_pass=-1)
    combined_run['metadata'][MI_DP]  = combine_metadata_element(data, MI_DP,  mixed_pass=-1)
    combined_run['metadata'][MI_CDP] = combine_metadata_element(data, MI_CDP, mixed_pass=-1)
    combined_run['metadata'][MI_OBS] = combine_metadata_element(data, MI_OBS, mixed_error=-1)

    # new UIDs from this
    combined_run['uid'] = run_uid_from_metadata(combined_run['metadata'])
    combined_run['uid_ignoring_fields'] = run_uid_from_metadata_ignoring_fields(combined_run['metadata'], fields_to_combine)

    # combine the observations as well
    observables = {observable for observable in run['observations']}
    for observable in observables:
        combined_run['observations'][observable] = combine_observation(observable, data)
    return combined_run

def parse_einsim_output_file(infilename):
    dirty_data = []
    ecc_code = None
    malformed_line = -1
    file_size_bytes = os.path.getsize(infilename)    
    with open_sim_file(infilename) as f:
        line_num = 0
        try:
            line = f.readline()
            while line:
                line_num += 1
                try:
                    parsed_line = parse_einsim_output_line(line, [], False)
                    if parsed_line:
                        dirty_data.append(parsed_line)
                    elif line.startswith('[ECC]'):
                        ecc_code = json.loads(line[6:])
                except KeyboardInterrupt:
                    raise
                except:
                    if malformed_line == -1:
                        print("[WARNING] malformed line in file " + infilename + " at line: " + str(line_num) + ' - ', sys.exc_info()[0])
                        print("    bad line:", line)
                        traceback.print_exc()
                        malformed_line = line_num
                    pass

                if line_num == 1 or (line_num % 1000) == 0:
                    sys.stdout.write('\t[' + infilename + "] {0:.2f}".format(100.0 * f.tell() / float(file_size_bytes)) + '% complete  \r')
                line = f.readline()
        except:
            print("Unexpected error in file " + infilename + " at line " + str(line_num) + ':', sys.exc_info()[0])
            raise

    return dirty_data, ecc_code

def run_uid_from_metadata_ignoring_fields(metadata, fields_to_combine):
    uid, nw, bl, bcl, ps, em, cd, dp, cdp, obs = metadata
    run_uid = '[DATA]'
    run_uid += ' uid:' + str(uid)       if 'uid'  not in fields_to_combine else ' uid:-1'
    run_uid += ' nw:' + str(nw)         if 'nw'   not in fields_to_combine else ' nw:-1'
    run_uid += ' bl:' + str(bl)         if 'bl'   not in fields_to_combine else ' bl:-1'
    run_uid += ' bcl:' + str(bcl)       if 'bcl'  not in fields_to_combine else ' bcl:-1'
    run_uid += ' ps:' + str(ps)         if 'ps'   not in fields_to_combine else ' ps:-1'
    run_uid += ' em:' + str(em)         if 'em'   not in fields_to_combine else ' em:-1'
    run_uid += ' cd:' + str(cd)         if 'cd'   not in fields_to_combine else ' cd:-1'
    run_uid += ' dp:' + str(dp)         if 'dp'   not in fields_to_combine else ' dp:-1'
    run_uid += ' cdp:' + str(cdp)       if 'cdp'  not in fields_to_combine else ' cdp:-1'
    run_uid += ' obs:' + str(obs)       if 'obs'  not in fields_to_combine else ' obs:-1'
    return run_uid

def run_uid_from_metadata(metadata):
    return run_uid_from_metadata_ignoring_fields(metadata, [])

# uses the UID, which contains the most information relative to the uid_ignoring_fields
def get_output_string(raw_data):
    out_str = ""
    assert raw_data['observations'], "Must have at least one observation in the data"
    for observable in raw_data['observations']:
        out_str += raw_data['uid']
        out_str += ' obs:' + observable
        if observable == OBS_N_ERRORS_PER_BURST:
            out_str += ' [ '
            out_str += ' '.join([':'.join(map(str, [pair[0], pair[1][0], pair[1][1]])) for pair in val['raw_data'].items()])
            out_str += ' ]\n'
        elif observable == OBS_PER_BIT_ERROR_COUNT:
            out_str += ' [ '
            for i in raw_data['observations'][observable]['hist_databurst']: out_str += str(i) + ' '
            out_str += ': '
            for i in raw_data['observations'][observable]['hist_codeburst']: out_str += str(i) + ' '
            out_str += ']\n'
    return out_str
    
# callback(run_uid, run_data, *callback_args) is called for every line parsed
# assumes that the input has already been aggregated! if not, repeat models will be treated as independent models
def parse_files_incremental(infilenames, is_experimental, fields_to_combine, callback, *callback_args):
    experiment_id = 0
    ecc_codes = {}

    # extract data from all files and aggregrate it
    print("parsing %d files" % len(infilenames))
    for infilename in infilenames:

        # for each line in the file, invoke the callback
        malformed_line = -1
        file_size_bytes = os.path.getsize(infilename)
        ecc_code = None
        with open_sim_file(infilename) as f:
            line_num = 0

            line = f.readline()
            while line:
                line_num += 1
                try:
                    parsed_data = parse_einsim_output_line(line, fields_to_combine, is_experimental)
                except KeyboardInterrupt:
                    raise
                except:
                    parsed_data = None
                    if malformed_line == -1:
                        print("[WARNING] malformed line in file " + infilename + " at line: " + str(line_num) + ' - ', sys.exc_info()[0])
                        print("    bad line:", line)
                        traceback.print_exc()
                        malformed_line = line_num
                    pass

                if parsed_data:
                    raw_data = parsed_data
                    callback(raw_data, ecc_code, *callback_args)
                elif line.startswith('[ECC]'):
                    ecc_codes[infilename] = json.loads(line[6:])

                if line_num == 1 or (line_num % 577) == 0:
                    sys.stdout.write('\t[' + infilename + "] {0:.2f}".format(100.0 * f.tell() / float(file_size_bytes)) + '% complete  \r')
                line = f.readline()
        sys.stdout.write('\n')
    return ecc_codes
