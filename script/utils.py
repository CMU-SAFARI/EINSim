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
import traceback
try:
    import lzma
except ImportError:
    print("[WARNING] unable to import lzma module required to read lzma/xz files. Use Python3.3+ or backports lzma for this functionality")
import gzip
import bz2

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
def parse_einsim_output_line(line):
    if len(line) > 6 and line[:6] == "[DATA]" in line:
        meta, raw_data = line[7:].strip().split('[')
        meta_split = meta.split(' ')
        ecc_category = meta_split[0][:-1]
        if ecc_category == 'REP':
            p = int(meta_split[1][2:])
            t = int(meta_split[2][2:])
            k = int(meta_split[3][2:])
            n = int(meta_split[4][2:])
            m = -1
            rber = float(meta_split[5][5:])
            bl = int(meta_split[6][3:])
            bcl = int(meta_split[7][4:])
            ps = int(meta_split[8][3:])
            ed = meta_split[9][3:]
            cd = meta_split[10][3:]
            dp = meta_split[11][3:]
        elif ecc_category == 'HSC':
            p = int(meta_split[1][2:])
            t = int(meta_split[2][2:])
            k = int(meta_split[3][2:])
            n = int(meta_split[4][2:])
            m = -1
            rber = float(meta_split[5][5:])
            bl = int(meta_split[6][3:])
            bcl = int(meta_split[7][4:])
            ps = int(meta_split[8][3:])
            ed = meta_split[9][3:]
            cd = meta_split[10][3:]
            dp = meta_split[11][3:]
        elif ecc_category == 'BCH':
            p = int(meta_split[1][2:])
            t = int(meta_split[2][2:])
            k = int(meta_split[3][2:])
            n = int(meta_split[4][2:])
            m = int(meta_split[5][2:])
            rber = float(meta_split[6][5:])
            bl = int(meta_split[7][3:])
            bcl = int(meta_split[8][4:])
            ps = int(meta_split[9][3:])
            ed = meta_split[10][3:]
            cd = meta_split[11][3:]
            dp = meta_split[12][3:]
        elif ecc_category == 'UNK':
            p = int(meta_split[1][2:])
            t = int(meta_split[2][2:])
            k = int(meta_split[3][2:])
            n = int(meta_split[4][2:])
            # only BCH has the extra m param
            if meta_split[5][0] == 'm':
                m = int(meta_split[5][2:])
                off = 1
            else:
                m = -1
                off = 0
            rber = float(meta_split[5 + off][5:])
            bl = int(meta_split[6 + off][3:])
            bcl = int(meta_split[7 + off][4:])
            ps = int(meta_split[8 + off][3:])
            ed = meta_split[9 + off][3:]
            cd = meta_split[10 + off][3:]
            dp = meta_split[11 + off][3:]
        else:
            print("[WARNING] unhandled ECC category \"" + ecc_category + "\" at line: " + str(line_num) + " in file " + infilename)
            traceback.print_exc()
            return None
        
        metadata = [ecc_category, p, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp]
        data = [map(int, entry.split(':')) for entry in raw_data.split(' ') if ':' in entry]

        # print(ecc_type, permutation, nerrs_correctable, n_data_bits, n_code_bits, dp, error_model)
        # print(hist_data)
        return [metadata, data]
    else:
        return None

def parse_einsim_output_file(infilename):
    dirty_data = []
    malformed_line = -1
    file_size_bytes = os.path.getsize(infilename)    
    with open_sim_file(infilename) as f:
        line_num = 0
        try:
            line = f.readline()
            while line:
                line_num += 1
                try:
                    parsed_line = parse_einsim_output_line(line)
                    if parsed_line:
                        dirty_data.append(parsed_line)
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

    return dirty_data

def run_uid_from_metadata(metadata):
    ecc_category, p, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp = metadata

    run_uid = '[DATA]'
    if ecc_category == 'REP':
        run_uid += ' REP:'
        run_uid += ' p:' + '-1'
        run_uid += ' t:' + str(t)
        run_uid += ' k:' + str(k)
        run_uid += ' n:' + str(n)
        run_uid += ' rber:' + str(rber)
        run_uid += ' bl:' + str(bl)
        run_uid += ' bcl:' + str(bcl)
        run_uid += ' ps:' + str(ps)
        run_uid += ' ed:' + str(ed)
        run_uid += ' cd:' + str(cd)
        run_uid += ' dp:' + str(dp)
    elif ecc_category == 'HSC':
        run_uid += ' HSC:'
        run_uid += ' p:' + '-1'
        run_uid += ' t:' + str(t)
        run_uid += ' k:' + str(k)
        run_uid += ' n:' + str(n)
        run_uid += ' rber:' + str(rber)
        run_uid += ' bl:' + str(bl)
        run_uid += ' bcl:' + str(bcl)
        run_uid += ' ps:' + str(ps)
        run_uid += ' ed:' + str(ed)
        run_uid += ' cd:' + str(cd)
        run_uid += ' dp:' + str(dp)
    elif ecc_category == 'BCH':
        run_uid += ' BCH:'
        run_uid += ' p:' + '-1'
        run_uid += ' t:' + str(t)
        run_uid += ' k:' + str(k)
        run_uid += ' n:' + str(n)
        run_uid += ' m:' + str(m)
        run_uid += ' rber:' + str(rber)
        run_uid += ' bl:' + str(bl)
        run_uid += ' bcl:' + str(bcl)
        run_uid += ' ps:' + str(ps)
        run_uid += ' ed:' + str(ed)
        run_uid += ' cd:' + str(cd)
        run_uid += ' dp:' + str(dp)
    elif ecc_category == 'UNK':
        run_uid += ' UNK:'
        run_uid += ' bl:' + str(bl)
    else:
        print("[WARNING] unhandled ECC category \"" + ecc_category)

    return run_uid

def parse_all_files(infilenames, experimental=False):
    raw_data = {}
    experiment_id = 0

    # extract data from all files and aggregrate it
    print("parsing %d files" % len(infilenames))
    for filename in infilenames:
        dirty_data = parse_einsim_output_file(filename)
        # print(dirty_data)
        for metadata, data in dirty_data:
            ecc_category, p, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp = metadata
            
            # ensure that experiments get unique UIDs - we do not aggregate them
            run_uid = run_uid_from_metadata(metadata)
            if experimental:
                run_uid += '_experiment_' + str(experiment_id)
                experiment_id += 1

            # omit the permutation - we're aggregating over it!
            all_metadata = [ecc_category, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp]
            if run_uid not in raw_data:
                raw_data[run_uid] = {'metadata' : all_metadata, 'data' : [], 'p' : []}
            raw_data[run_uid]['data'].append(data)
            raw_data[run_uid]['p'].append(p)

    # print("Done parsing all files. Got:", len(raw_data), "files worth")

    # combine histograms
    for run_uid, components in raw_data.items():
        ecc_category, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp = components['metadata']
        
        raw_histogram = {}
        re_histogram = {}
        ue_histogram = {}
        for run in components['data']:
            for n_errs, re_count, ue_count in run:
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
        raw_data[run_uid] = {'metadata' : components['metadata']
                            , 'p' : components['p']
                            , 'raw_data' : raw_histogram
                            , 'rber' : rber_calculated
                            , 'uber' : uber_calculated
                            , 're_data' : re_histogram
                            , 'ue_data': ue_histogram}

    return raw_data

def get_minimized_output_string(raw_data):
    out_str = ""
    for key, val in raw_data.items(): 
        out_str += key
        out_str += ' [ '
        out_str += ' '.join([':'.join(map(str, [pair[0], pair[1][0], pair[1][1]])) for pair in val['raw_data'].items()])
        out_str += ' ]\n'
    return out_str
    
# callback(run_uid, run_data, *callback_args) is called for every line parsed
# assumes that the input has already been aggregated! if not, repeat models will be treated as independent models
def parse_files_incremental(infilenames, is_experimental, callback, *callback_args):
    experiment_id = 0

    # extract data from all files and aggregrate it
    print("parsing %d files" % len(infilenames))
    for infilename in infilenames:

        # for each line in the file, invoke the callback
        malformed_line = -1
        file_size_bytes = os.path.getsize(infilename)
        with open_sim_file(infilename) as f:
            line_num = 0

            line = f.readline()
            while line:
                line_num += 1
                try:
                    parsed_line = parse_einsim_output_line(line)
                except KeyboardInterrupt:
                    raise
                except:
                    parsed_line = None
                    if malformed_line == -1:
                        print("[WARNING] malformed line in file " + infilename + " at line: " + str(line_num) + ' - ', sys.exc_info()[0])
                        print("    bad line:", line)
                        traceback.print_exc()
                        malformed_line = line_num
                    pass

                if parsed_line:
                    metadata, data = parsed_line
                    ecc_category, p, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp = metadata
                    
                    # ensure that experiments get unique UIDs - we do not aggregate them
                    run_uid = run_uid_from_metadata(metadata)
                    if is_experimental:
                        run_uid += '_experiment_' + str(experiment_id)
                        experiment_id += 1

                    # omit the permutation - we're aggregating over it!
                    metadata_no_p = [ecc_category, t, k, n, m, rber, bl, bcl, ps, ed, cd, dp]

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
                    raw_data = {'metadata' : metadata_no_p
                               , 'p' : [p]
                               , 'raw_data' : raw_histogram
                               , 'rber' : rber_calculated
                               , 'uber' : uber_calculated
                               , 're_data' : re_histogram
                               , 'ue_data': ue_histogram}
                    callback(run_uid, raw_data, *callback_args)

                if line_num == 1 or (line_num % 1000) == 0:
                    sys.stdout.write('\t[' + infilename + "] {0:.2f}".format(100.0 * f.tell() / float(file_size_bytes)) + '% complete  \r')
                line = f.readline()
        sys.stdout.write('\n')
