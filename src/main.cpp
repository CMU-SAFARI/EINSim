/**
 * @file main.cpp
 *
 * @brief Main file that determines what to do based on CLI arguments
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <cstdio>
#include <cinttypes>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <stdarg.h>

/* libraries */
#include "libtp/libtp.h"
#define CXXOPTS_VECTOR_DELIMITER ';' /**< hideous hack to allow backwards compatability */
#include "cxxopts/cxxopts.h"

/* project includes */
#include "word_generator.h"
#include "observable.h"
#include "error_model.h"
#include "ecc_code.h"
#include "simulate.h"
#include "debug.h"
#include "supporting_routines.h"
#include "codes/repetition_code.h"
#include "codes/bch_code.h"
#include "codes/hamming_code.h"

/**
 * @brief EINSim entry point that parses CLI arguments and begins simulation
 * 
 * @param argc number of command line arguments
 * @param argv command line arguments
 * 
 * @return int return code
 */
int main(int argc, char ** argv) // canNOT const argv due to cxxopts requirements
{
    // save the comand line as a string
    std::stringstream command_line;
    for(int i = 0; i < argc; i++)
        command_line << ' ' << argv[i];

    // parse the CLI options
    cxxopts::Options option_parser("einsim", "Probabilistic ECC simulator");
    option_parser.add_options("Common")
        ("m, mode", "Choose mode: (t)est, (d)ebug, or (s)imulate", cxxopts::value< std::string >())
        ("n, nwords", "# words to simulate", cxxopts::value< uint64_t >())
        ("x, max_words", "maximum # words to simulate per job", cxxopts::value< uint64_t >())
        ("t, nthreads", "# worker threads", cxxopts::value< int >())
        ("v, verbose", "Print non-essential messages")
        ("f, file", "Output file name", cxxopts::value< std::string >())
        ("h, help", "Show help")
        ;
    option_parser.add_options("Simulation")
        ("b, burst_length_bits", "Burst lengths to simulate (#data bits)", cxxopts::value< std::vector< int > >())
        ("w, word_to_burst_mapping", std::string("Mapping from individual ECC words to a DRAM burst {") + einsim::get_all_possible_word_to_burst_mappings() + "}", cxxopts::value< std::vector< std::string > >())
        ("c, true_anti_cell_distributions", std::string("True- and anti-cell distribution {") + einsim::get_all_possible_true_anti_cell_distributions() + "}", cxxopts::value< std::vector< std::string > >())
        ("d, data_patterns", std::string("Data pattern to simulate {") + einsim::get_all_possible_data_patterns() + "} OR custom (0b, 0o, 0x)", cxxopts::value< std::vector< std::string > >())
        ("e, error_models", std::string("Error model to use {") + einsim::get_all_possible_error_models() + "} (supply N comma-separated model specification tuples name0,p0,..,pn,name1,p0,... for N bits or just one for all bits) OR filename for error model JSON configuration file", cxxopts::value< std::vector< std::string > >())
        // ("i, fault_injection_mask", "fault-injection bit-mask of burst codeword length (0b..., 0o..., 0x...)", cxxopts::value< std::vector< std::string > >())
        ("o, observables", std::string("observables to measure from the data {") + einsim::get_all_possible_observables() + "}", cxxopts::value< std::vector< std::string > >())
        // ("r, rber", "RBER(s) to simulate (0 <= r <= 1)", cxxopts::value< std::vector< float > >())
        ("s, ecc_scheme", std::string("ECC scheme(s) to simulate {") + einsim::get_all_possible_ecc_schemes() +  "} OR filename for ECC scheme JSON configuration file", cxxopts::value< std::vector< std::string > >())
        ("k, data_bits", "number of ECC data bits to simulate (k >= 1)", cxxopts::value< std::vector< int > >())
        ("p, permutations", "permutations to compute (p => 0) [int || int-int]", cxxopts::value< std::vector< std::string > >())
        ("y, dry_run", "exit after printing configuration")
        ;
    option_parser.add_options("Test")
        ("T, test_mode", std::string("Test mode(s) to run {") + einsim::get_all_possible_test_modes() +  "}", cxxopts::value< std::vector< std::string > >())
        ;
    bool needs_help = (argc == 1);
    auto options = option_parser.parse(argc, argv);
    if(needs_help || options.count("help"))
    {
        std::cout << option_parser.help({"", "Common", "Test", "Simulation"}) << std::endl;
        return 0;
    }

    // set g_verbosity
    g_verbosity = options.count("verbose");

    // prepare Eigen library to use multithreading
    Eigen::initParallel();

    // check if output file is requested
    std::string output_filename;
    if(options.count("file"))
    {
        output_filename = options["file"].as< std::string >();
        printf("Redirecting output to file: \"%s\"\n", output_filename.c_str());

        // check if the file exists
        FILE *temp_file;
        if((temp_file = fopen(output_filename.c_str(), "r")))
        {
            fclose(temp_file);
            printf("[WARN] output file \"%s\" already exists!\n", output_filename.c_str());
            printf("[INFO] type 'y' to overwrite, any other key to exit\n");
            char key = std::cin.get();
            if(key == 'y')
            {
                if(std::remove(output_filename.c_str()) != 0)
                {
                    std::perror((std::string("Error deleting file ") + output_filename).c_str());
                    exit(-1);
                }
            }
            else
                exit(-1);
        }

        if(!(g_output_file = fopen(output_filename.c_str(), "wb")))
        {
            printf("[ERROR] output file \"%s\" could not be opened for writing!\n", output_filename.c_str());
            exit(-1);
        }
    }
    else
    {
        g_output_file = stdout;
        printf("[WARNING] No output file specified - using only stdout\n");
    }

    // echo the command line
    printf_both("[INFO] executable command:%s\n", command_line.str().c_str());

    uint64_t n_bursts_to_simulate = 100; 
    if(options.count("nwords"))
        n_bursts_to_simulate = options["nwords"].as< uint64_t >();

    int n_worker_threads = 1;
    if(options.count("nthreads"))
    {
        n_worker_threads = options["nthreads"].as< int >();
        printf_both("[INFO] using %d threads\n", n_worker_threads);
    }
    else
        printf_both("[WARNING] no thread count specified- using %d threads\n", n_worker_threads);

    // handle mode-specific details
    if(options.count("mode") != 1)
    {
        printf("[ERROR] Must choose exactly one mode\n");
        std::cout << option_parser.help({"", "Common", "Test", "Simulation"}) << std::endl;
        return -1;
    }

    std::string mode = options["mode"].as< std::string >();
    if(mode == "t")
    {
        const char *status = "[INFO] Configuring for Test Mode";
        if(g_verbosity > 0)
            printf_both("%s\n", status);
        else
            fprintf(g_output_file, "%s\n", status);

        // get the the test mode option
        if(options.count("test_mode") == 0)
        {
            printf("[ERROR] Must provide test modes to simulate\n");
            std::cout << option_parser.help({"", "Test"}) << std::endl;
            return -1;
        }
        const std::vector< std::string > test_mode_strs = options["test_mode"].as< std::vector< std::string > >();
        std::set< enum einsim::test_mode > test_modes;
        for(const std::string &test_mode_str : test_mode_strs)
        {
            enum einsim::test_mode tm = einsim::str_to_enum_test_mode(test_mode_str);
            if(tm == einsim::TM_UNKNOWN)
            {
                printf("[ERROR] Invalid test_mode: %s\n", test_mode_str.c_str());
                std::cout << option_parser.help({"", "Test"}) << std::endl;
                return -1;            
            }
            test_modes.insert(tm);
        }

        // run the test for all ECC schemes
        for(const enum einsim::test_mode tm : test_modes)
        {
            einsim::test_ecc(einsim::hamming::submit_tests, tm, n_worker_threads);
            einsim::test_ecc(einsim::bch::submit_tests, tm, n_worker_threads);
            einsim::test_ecc(einsim::repetition::submit_tests, tm, n_worker_threads);
        }
    }
    else if(mode == "d")
    {
        const char *status = "[INFO] Configuring for Debug Mode";
        if(g_verbosity > 0)
            printf_both("%s\n", status);
        else
            fprintf(g_output_file, "%s\n", status);

        debug_example(n_worker_threads, n_bursts_to_simulate);
    }
    else if(mode == "s")
    {
        const char *status = "[INFO] Configuring for Simulation Mode";
        if(g_verbosity > 0)
            printf_both("%s\n", status);
        else
            fprintf(g_output_file, "%s\n", status);

/**
 * @brief generic error message to print alongside the usage command
 * 
 * @param msg message to print
 */
#define SIMULATION_OPTION_REQUIRED(msg) do {                                \
        printf(msg); \
        std::cout << option_parser.help({"", "Simulation"}) << std::endl;   \
        return -1;                                                          \
    } while(0)

        /************************************************************************
         * determine the simulation parameters from the CLI
         ************************************************************************/
        // burst-length
        std::set< int > burst_lengths;
        if(options.count("burst_length_bits") == 0)
            SIMULATION_OPTION_REQUIRED("[ERROR] must provide at least one burst length to simulate\n");
        else
        {
            const std::vector< int > &burst_lengths_raw = options["burst_length_bits"].as< std::vector< int > >();
            burst_lengths.insert(burst_lengths_raw.begin(), burst_lengths_raw.end());
        }

        // ecc word to burst mapping
        std::set< enum einsim::word_to_burst_mapping > w2b_mappings;
        if(options.count("word_to_burst_mapping") == 0)
            w2b_mappings = std::set< enum einsim::word_to_burst_mapping >({einsim::W2BM_BLOCKS});
        else
        {
            const std::vector< std::string > &w2b_mappings_vec = options["word_to_burst_mapping"].as< std::vector< std::string > >();
            for(const std::string &w2b_mapping_str : w2b_mappings_vec)
            {
                enum einsim::word_to_burst_mapping w2bm = einsim::str_to_enum_word_to_burst_mapping(w2b_mapping_str);
                if(w2bm == einsim::W2BM_UNKNOWN)
                {
                    printf("[ERROR] Invalid word-to-burst mapping: %s\n", w2b_mapping_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;
                }
                w2b_mappings.insert(w2bm);
            }
        }

        // data-pattern
        std::vector< enum einsim::data_pattern > data_patterns;
        std::vector< Eigen::Matrix< ET, Eigen::Dynamic, 1 > > custom_dps;
        if(options.count("data_patterns") == 0)
            SIMULATION_OPTION_REQUIRED("[ERROR] must provide at least one data pattern to simulate\n");
        else
        {
            const std::vector< std::string > dp_list = options["data_patterns"].as< std::vector< std::string > >();
            for(const std::string &data_pattern_str : dp_list)
            {
                enum einsim::data_pattern dp = einsim::str_to_enum_data_pattern(data_pattern_str);
                if(dp == einsim::DP_CUSTOM)
                    custom_dps.push_back(einsim::custom_dp_to_vector(data_pattern_str));
                else if(dp == einsim::DP_UNKNOWN)
                {
                    printf("[ERROR] Invalid data pattern: %s\n", data_pattern_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;            
                }
                data_patterns.push_back(dp);
            }
        }

        // true-/anti-cell distribution 
        std::set< enum einsim::true_anti_cell_distribution > true_anti_cell_distributions;
        if(options.count("true_anti_cell_distributions") == 0)
            true_anti_cell_distributions = std::set< enum einsim::true_anti_cell_distribution >({einsim::CD_ALL_TRUE_OR_ALL_ANTI});
        else
        {
            const std::vector< std::string > cd_list = options["true_anti_cell_distributions"].as< std::vector< std::string > >();
            for(const std::string &true_anti_cell_distribution_str : cd_list)
            {
                enum einsim::true_anti_cell_distribution cd = einsim::str_to_enum_true_anti_cell_distribution(true_anti_cell_distribution_str);
                if(cd == einsim::CD_UNKNOWN)
                {
                    printf("[ERROR] Invalid true-/anti-cell distribution: %s\n", true_anti_cell_distribution_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;            
                }
                true_anti_cell_distributions.insert(cd);
            }
        }

        // observables - we currently only implement one, so we will always use it
        std::set< enum einsim::observable > observables;
        if(options.count("observables") == 0)
            SIMULATION_OPTION_REQUIRED("[ERROR] must provide at least one observable to simulate\n");
        else
        {
            const std::vector< std::string > cd_list = options["observables"].as< std::vector< std::string > >();
            for(const std::string &observable_str : cd_list)
            {
                enum einsim::observable cd = einsim::str_to_enum_observable(observable_str);
                if(cd == einsim::OBS_UNKNOWN)
                {
                    printf("[ERROR] Invalid observable: %s\n", observable_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;            
                }
                observables.insert(cd);
            }
        }

        // fault injection mask
        // std::vector< Eigen::Matrix< ET, Eigen::Dynamic, 1 > > fault_injection_mask_list;
        // if(options.count("fault_injection_mask") != 0)
        // {
        //     const std::vector< std::string > fault_injection_masks = options["fault_injection_mask"].as< std::vector< std::string > >();
        //     for(const std::string &fault_injection_mask_str : fault_injection_masks)
        //         fault_injection_mask_list.push_back(einsim::custom_fault_injection_mask_to_vector(fault_injection_mask_str));
        // }

        // set the number of words per simulation job
        uint64_t n_bursts_per_job = 10000;
        if(options.count("max_words"))
            n_bursts_per_job = options["max_words"].as< uint64_t >();

        /*************************************
         * generate some ECC schemes
         *************************************/ 
        // n data bits
        std::set< int > n_ecc_data_bits_parameterized; // = {4, 8, 16, 32, 64, 128, 256};
        if(options.count("data_bits") != 0)
        {
            const std::vector< int > k_list = options["data_bits"].as< std::vector< int > >();
            n_ecc_data_bits_parameterized.insert(k_list.begin(), k_list.end());
        }
  
        // permutations
        std::set< int > permutations_parameterized; // = {0, 1, 2...}
        if(options.count("permutations") != 0)
        {
            const std::vector< std::string > p_list = options["permutations"].as< std::vector< std::string > >();
            for(const std::string& entry : p_list)
            {
                std::size_t found = entry.find("-");
                if(found != std::string::npos)
                {
                    int64_t range_start = std::stoll(entry.substr(0, found));
                    int64_t range_end = std::stoll(entry.substr(found + 1));
                    // std::cout << "range:" << range_start << " to " << range_end << '\n';
                    assert(range_start <= range_end && "range list must be in order");
                    for(int64_t i = range_start; i <= range_end; i++)
                        permutations_parameterized.insert(i);
                }
                else
                    permutations_parameterized.insert(std::stoll(entry));
            }
        }
        
        // determine how many ECC schemes we need to build - either CLI-based (parameterized) or a config file
        std::vector< enum einsim::ecc_scheme > ecc_schemes_parameterized;
        std::vector< std::string > ecc_scheme_cfg_filenames;
        if(options.count("ecc_scheme") == 0)
            SIMULATION_OPTION_REQUIRED("ecc_scheme");
        else
        {
            const std::vector< std::string > &scheme_list = options["ecc_scheme"].as< std::vector< std::string > >();
            for(const std::string &ecc_scheme_str : scheme_list)
            {
                // try to open a file
                std::ifstream f(ecc_scheme_str, std::ifstream::in);
                if(f.good())
                {
                    ecc_scheme_cfg_filenames.push_back(ecc_scheme_str);
                    f.close();
                }
                else
                {
                    // instantiate the different ECC types for the data bit sizes and ECC schemes requested
                    enum einsim::ecc_scheme ecc_scheme = einsim::str_to_enum_ecc_scheme(ecc_scheme_str);
                    if(ecc_scheme == einsim::ES_UNKNOWN)
                    {
                        printf("[ERROR] unknown/invalid ECC scheme requested: %s\n", ecc_scheme_str.c_str());
                        std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                        return -1;            
                    }
                    else
                        ecc_schemes_parameterized.push_back(ecc_scheme);
                }
            }
        }

        // build the ECC schemes
        std::vector< einsim::ecc_code * > ecc_schemes;
        std::set< enum einsim::ecc_scheme > ecc_schemes_cfg_files;
        std::set< int > n_ecc_data_bits_cfg_files;
        std::set< int > permutations_cfg_files;
        if(ecc_scheme_cfg_filenames.size())
        {
            std::cout << "[INFO] building " << ecc_scheme_cfg_filenames.size() << " ECC schemes from configuration files" << std::endl;
            for(const std::string &cfg_file_name : ecc_scheme_cfg_filenames)
            {
                einsim::ecc_code *code = einsim::build_ecc_code(cfg_file_name);
                if(code == NULL)
                {
                    std::cout << "[ERROR] unable to build ECC code for configuration file: " << cfg_file_name << std::endl;
                    return -1;
                }
                ecc_schemes.push_back(code);
                ecc_schemes_cfg_files.insert(code->get_scheme());
                n_ecc_data_bits_cfg_files.insert(code->get_n_data_bits());
                permutations_cfg_files.insert(code->get_permutation());
            }
        }
        if(ecc_schemes_parameterized.size())
        {
            std::cout << "[INFO] building " << ecc_schemes_parameterized.size() << " ECC schemes over " << permutations_parameterized.size() * n_ecc_data_bits_parameterized.size() 
                << " configurations each (" << ecc_schemes_parameterized.size() * permutations_parameterized.size() * n_ecc_data_bits_parameterized.size() << ") schemes total" << std::endl;
            if(permutations_parameterized.size() == 0)
                SIMULATION_OPTION_REQUIRED("[ERROR] must specify at least one permutation when specifying an ECC scheme using ECC code parameters");
            for(int ecc_permutation : permutations_parameterized)
            {
                if(n_ecc_data_bits_parameterized.size() == 0)
                    SIMULATION_OPTION_REQUIRED("[ERROR] must specify at least one data bit size when specifying an ECC scheme using ECC code parameters");
                for(int n_data_bits : n_ecc_data_bits_parameterized)
                    for(einsim::ecc_scheme scheme : ecc_schemes_parameterized)
                    {
                        einsim::ecc_code *code = einsim::build_ecc_code(scheme, n_data_bits, ecc_permutation);
                        if(code == NULL)
                        {
                            std::cout << "[ERROR] unable to build ECC code for configuration p: " << ecc_permutation << " k: " << n_data_bits 
                                << " s: " << enum_to_str_ecc_scheme(scheme) << std::endl;
                            return -1;
                        }
                        ecc_schemes.push_back(code);
                    }
            }
        }

        // build the error models
        // each argument is one model (or a range)
        // supplying N parameters for a model is interpreted as one-per-bit, and one parameter is N bits
        // e.g., -e UNIFORM_RANDOM,0.1,0.2,0.3,STUCK_AT,0,1 is TWO bits with the cartesian product of 3x2 models
        std::vector< std::vector< einsim::error_model_descriptor * > > error_models;
        if(options.count("error_models") == 0)
            SIMULATION_OPTION_REQUIRED("error_models");
        else
        {
            const std::vector< std::string > em_list = options["error_models"].as< std::vector< std::string > >();
            for(const std::string &em_list_entry : em_list)
            {
                // try to open a file
                std::ifstream f(em_list_entry, std::ifstream::in);
                if(f.good())
                {
                    const std::string &cfg_file_name = em_list_entry;
                    f.close();
                    std::cout << "[INFO] building error model from configuration file: " << cfg_file_name << std::endl;
                    std::vector< std::vector< einsim::error_model_descriptor * > > emd_vec_vec = einsim::error_model_descriptors_from_json(cfg_file_name);
                    for(const std::vector< einsim::error_model_descriptor * > &emd_vec : emd_vec_vec)
                        error_models.push_back(emd_vec);
                }
                else
                {
                    // must be a comma separated string, separate the commas
                    std::stringstream ss(em_list_entry);
                    std::vector< std::string > em_list_entry_split;
                    while(ss.good())
                    {
                        std::string substr;
                        std::getline(ss, substr, ',');
                        em_list_entry_split.push_back(substr);
                    }

                    // build a vector of the different bit models!
                    std::vector< std::vector< einsim::error_model_descriptor * > > emds_per_bit;
                    for(size_t i = 0; i < em_list_entry_split.size();)
                    {
                        const std::string &em_list_token = em_list_entry_split.at(i);
                        // std::cout << "[DEBUG] " << i << " " << em_list_token << std::endl;
                        enum einsim::error_model em = einsim::str_to_enum_error_model(em_list_token);
                        int n_model_params = einsim::error_model_descriptor::get_n_model_params(em);
                        if(n_model_params == -1 || (i + n_model_params >= em_list_entry_split.size()))
                        {
                            std::cout << "[ERROR] incorrect number of model parameters given for error model " << em_list_token 
                                << " - expected: " << n_model_params 
                                << " got: " << (em_list_entry_split.size() - i - n_model_params)  << std::endl;
                            std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                            return -1;               
                        }
                        // std::cout << "[DEBUG] " << i << " " << enum_to_str_error_model(em) << " " << n_model_params << std::endl;

                        i++;
                        std::vector< einsim::error_model_descriptor * > emds_this_bit;
                        while(true)
                        {
                            std::vector< std::string > model_params(em_list_entry_split.begin() + i, em_list_entry_split.begin() + i + n_model_params);
                            i += n_model_params;
                            einsim::error_model_descriptor *emd = einsim::error_model_descriptor_from_params(em, model_params);
                            emds_this_bit.push_back(emd);

                            if(i == em_list_entry_split.size())
                                break;
                            else
                            {
                                const std::string &next_token = em_list_entry_split.at(i);
                                enum einsim::error_model em = einsim::str_to_enum_error_model(next_token);
                                if(em != einsim::EM_UNKNOWN)
                                    break;
                            }
                        }
                        emds_per_bit.push_back(emds_this_bit);
                    }

                    // construct the cross-product space of these models
                    einsim::construct_cartesian_product_of_per_bit_error_models(error_models, emds_per_bit);
                }
            }
        }

        // check UIDs and ensure that they are all unique
        std::set< uint64_t > ecc_scheme_uids;
        for(const einsim::ecc_code *ec : ecc_schemes)
            ecc_scheme_uids.insert(ec->get_uid());
        assert(ecc_scheme_uids.size() == ecc_schemes.size() && "UID hash collision detected!");
        
        /************************************************************************
         * launch the desired configuration
         ************************************************************************/
        uint64_t n_configs = burst_lengths.size() * w2b_mappings.size() * data_patterns.size() 
            * error_models.size() * true_anti_cell_distributions.size() * observables.size() 
            * ecc_schemes.size();
        std::cout << "[INFO] testing " << n_configs << " configurations subdivided into groups of " << n_bursts_per_job << " bursts per job:" << std::endl;
        std::cout << "[INFO]    " << burst_lengths.size()                << " burst_length_bits:            [ "; for(const auto &i : burst_lengths)                std::cout << i << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << w2b_mappings.size()                 << " word-to-burst mappings:       [ "; for(const auto &i : w2b_mappings)                 std::cout << einsim::enum_to_str_word_to_burst_mapping(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << data_patterns.size()                << " data_patterns:                [ "; for(const auto &i : data_patterns)                std::cout << einsim::enum_to_str_data_pattern(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << custom_dps.size()               << " custom_patterns:" << std::endl; 
        for(const auto &i : custom_dps) 
            std::cout << "[INFO]        [" << i.transpose()          << ']' << std::endl;
        std::cout << "[INFO]    " << error_models.size()          << " error_models:" << std::endl;
        for(const std::vector< einsim::error_model_descriptor * > &emd_vec : error_models)
        {
            std::cout << "[INFO]        ["; 
            for(const einsim::error_model_descriptor *emd : emd_vec) 
                std::cout << emd->to_str() << " "; 
            std::cout << "]" << std::endl;
        }
        std::cout << "[INFO]    " << true_anti_cell_distributions.size() << " true_anti_cell_distributions: [ "; for(const auto &i : true_anti_cell_distributions) std::cout << einsim::enum_to_str_true_anti_cell_distribution(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << observables.size()                  << " observables:                  [ "; for(const auto &i : observables)                  std::cout << einsim::enum_to_str_observable(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << ecc_schemes.size()                    << " ECC schemes:" << std::endl;
        std::cout << "[INFO]        generated from code parameters:" << std::endl;
        std::cout << "[INFO]            " << ecc_schemes_parameterized.size()     << " schemes:                      [ "; for(const auto &i : ecc_schemes_parameterized)     std::cout << einsim::enum_to_str_ecc_scheme(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]            " << n_ecc_data_bits_parameterized.size() << " data-words:                   [ "; for(const auto &i : n_ecc_data_bits_parameterized) std::cout << i << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]            " << permutations_parameterized.size()    << " permutations:                 [ "; print_ranges(permutations_parameterized);
        std::cout << "[INFO]        read from cfg files:" << std::endl;
        std::cout << "[INFO]            " << ecc_schemes_cfg_files.size()         << " schemes:                      [ "; for(const auto &i : ecc_schemes_cfg_files)         std::cout << einsim::enum_to_str_ecc_scheme(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]            " << n_ecc_data_bits_cfg_files.size()     << " data-words:                   [ "; for(const auto &i : n_ecc_data_bits_cfg_files)     std::cout << i << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]            " << permutations_cfg_files.size()        << " permutations:                 [ "; print_ranges(permutations_cfg_files);
        
        if(options.count("dry_run") != 0)
        {
            std::cout << "[INFO] dry run complete" << std::endl;
            exit(0);
        }

        // run the requested simulation configuration
        simulate
        (
              n_worker_threads
            , n_bursts_to_simulate
            , n_bursts_per_job
            , burst_lengths
            , w2b_mappings
            , data_patterns
            , custom_dps
            , error_models
            , true_anti_cell_distributions
            , observables
            , ecc_schemes
        );
    }
    else
    {
        printf("[ERROR] invalid mode!\n");
        std::cout << option_parser.help({"", "Common", "Test", "Simulation"}) << std::endl;
        return -1;
    }

    return 0;
}
