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
#include <stdarg.h>

/* libraries */
#include "libtp.h"
#include "cxxopts.h"

/* project includes */
#include "word_generator.h"
#include "observable.h"
#include "error_injector.h"
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
 * @return return code
 */
int main(int argc, const char **argv)
{
    cxxopts::Options option_parser("einsim", "Probabilistic ECC simulator");
    option_parser.add_options("Common")
        ("m, mode", "Choose mode: (t)est, (d)ebug, or (s)imulate", cxxopts::value< std::string >())
        ("n, nwords", "# words to simulate", cxxopts::value< int >())
        ("t, nthreads", "# worker threads", cxxopts::value< int >())
        ("v, verbose", "Print non-essential messages")
        ("o, output", "Output file name", cxxopts::value< std::string >())
        ("h, help", "Show help")
        ;
    option_parser.add_options("Simulation")
        ("b, burst_length_bits", "Burst length to simulate (bits) (b >= 0)", cxxopts::value< int >())
        ("d, data_patterns", std::string("Data pattern to simulate {") + einsim::get_all_possible_data_patterns() + "}", cxxopts::value< std::vector< std::string > >())
        ("e, error_distributions", std::string("Error distribution to model {") + einsim::get_all_possible_error_distributions() + "}", cxxopts::value< std::vector< std::string > >())
        ("c, true_anti_cell_distribution", std::string("True- and anti-cell distribution {") + einsim::get_all_possible_true_anti_cell_distributions() + "}", cxxopts::value< std::vector< std::string > >())
        ("r, rber", "RBER(s) to simulate (0 <= r <= 1)", cxxopts::value< std::vector< float > >())
        ("k, data_bits", "number of ECC data bits to simulate (k >= 1)", cxxopts::value< std::vector< int > >())
        ("s, ecc_scheme", std::string("ECC scheme(s) to simulate {") + einsim::get_all_possible_ecc_schemes() +  "}", cxxopts::value< std::vector< std::string > >())
        ("p, permutations", "permutations to compute (p => 0)", cxxopts::value< std::vector< int > >())
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

    // check if output file is requested
    std::string output_filename;
    if(options.count("output"))
    {
        output_filename = options["output"].as< std::string >();
        printf("Redirecting output to file: \"%s\"\n", output_filename.c_str());

        // check if the file exists
        FILE *temp_file;
        if((temp_file = fopen(output_filename.c_str(), "r")))
        {
            fclose(temp_file);
            printf("[ERROR] output file \"%s\" already exists!\n", output_filename.c_str());
            exit(-1);
        }
        else
        {
            if(!(g_output_file = fopen(output_filename.c_str(), "wb")))
            {
                printf("[ERROR] output file \"%s\" could not be opened for writing!\n", output_filename.c_str());
                exit(-1);
            }
        }
    }
    else
    {
        g_output_file = stdout;
        printf("[WARNING] No output file specified - using only stdout\n");
    }

    int n_words_to_simulate = 100; 
    if(options.count("nwords"))
        n_words_to_simulate = options["nwords"].as< int >();

    int n_worker_threads = 1;
    if(options.count("nthreads"))
    {
        n_worker_threads = options["nthreads"].as< int >();
        printf_both("using %d threads\n", n_worker_threads);
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
        for(std::string test_mode_str : test_mode_strs)
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

        debug_example(n_worker_threads, n_words_to_simulate);
    }
    else if(mode == "s")
    {
        const char *status = "[INFO] Configuring for Simulation Mode";
        if(g_verbosity > 0)
            printf_both("%s\n", status);
        else
            fprintf(g_output_file, "%s\n", status);

        /************************************************************************
         * determine the simulation parameters from the CLI
         ************************************************************************/
        // burst-length
        if(options.count("burst_length_bits") == 0)
        {
            printf("[ERROR] Must provide exactly one burst length to simulate\n");
            std::cout << option_parser.help({"", "Simulation"}) << std::endl;
            return -1;
        }
        int bl = options["burst_length_bits"].as< int >();

        // data-pattern
        std::set< enum einsim::data_pattern > data_patterns;
        if(options.count("data_patterns") == 0)
        {
            data_patterns = std::set< enum einsim::data_pattern >({einsim::DP_RANDOM});
        }
        else
        {
            const std::vector< std::string > dp_list = options["data_patterns"].as< std::vector< std::string > >();
            for(std::string data_pattern_str : dp_list)
            {
                enum einsim::data_pattern dp = einsim::str_to_enum_data_pattern(data_pattern_str);
                if(dp == einsim::DP_UNKNOWN)
                {
                    printf("[ERROR] Invalid data pattern: %s\n", data_pattern_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;            
                }
                data_patterns.insert(dp);
            }
        }

        // error-distribution
        std::set< enum einsim::error_distribution > error_distributions;
        if(options.count("error_distributions") == 0)
        {
            error_distributions = std::set< enum einsim::error_distribution >({einsim::ED_UNIFORM_RANDOM});
        }
        else
        {
            const std::vector< std::string > ed_list = options["error_distributions"].as< std::vector< std::string > >();
            for(std::string error_distribution_str : ed_list)
            {
                enum einsim::error_distribution ed = einsim::str_to_enum_error_distribution(error_distribution_str);
                if(ed == einsim::ED_UNKNOWN)
                {
                    printf("[ERROR] Invalid error distribution: %s\n", error_distribution_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;            
                }
                error_distributions.insert(ed);
            }
        }

        // true-/anti-cell distribution 
        std::set< enum einsim::true_anti_cell_distribution > true_anti_cell_distributions;
        if(options.count("true_anti_cell_distributions") == 0)
        {
            true_anti_cell_distributions = std::set< enum einsim::true_anti_cell_distribution >({einsim::CD_ALL_TRUE_OR_ALL_ANTI});
        }
        else
        {
            const std::vector< std::string > cd_list = options["true_anti_cell_distributions"].as< std::vector< std::string > >();
            for(std::string true_anti_cell_distribution_str : cd_list)
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

        // RBER
        std::set< float > rbers;//  = {1e-6, 2.5e-6, 5e-6, 7.5e-6, 1e-5, 2.5e-5, 5e-5, 7.5e-5, 1e-4, 2.5e-4, 5e-4, 7.5e-4, 1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 2.5e-2, 5e-2, 7.5e-2, 1e-1};
        if(options.count("rber") == 0)
        {
            for(float rber = 1e-6; rber < 0.5; rber *= 1.01)
                rbers.insert(rber);
        }
        else
        {
            const std::vector< float > rber_list = options["rber"].as< std::vector< float > >();
            rbers.insert(rber_list.begin(), rber_list.end());
        }
        
        std::set< int > n_ecc_data_bits; // = {4, 8, 16, 32, 64, 128, 256};
        if(options.count("data_bits") == 0)
        {
            n_ecc_data_bits = std::set< int >({4, 8, 16, 32, 64, 128, 256});
        }
        else
        {
            const std::vector< int > k_list = options["data_bits"].as< std::vector< int > >();
            n_ecc_data_bits.insert(k_list.begin(), k_list.end());
        }
  
        // permutations
        std::set< int > permutations; // = {0, 1, 2...};
        if(options.count("permutations") == 0)
        {
            std::mt19937 rng;
            rng.seed(std::random_device()());
            std::uniform_int_distribution<std::mt19937::result_type> dist(0, -1u);
            permutations = std::set< int >({(int)dist(rng)});
        }
        else
        {
            const std::vector< int > p_list = options["permutations"].as< std::vector< int > >();
            permutations.insert(p_list.begin(), p_list.end());
        }
        
        // ECC schemes
        std::set< enum einsim::ecc_scheme > ecc_schemes; // = {einsim::REPETITION_T1, einsim::REPETITION_T2, einsim::REPETITION_T3, einsim::HAMMING_SEC, einsim::BCH_T1, einsim::BCH_T2, einsim::BCH_T3};
        if(options.count("ecc_scheme") == 0)
        {
            ecc_schemes = std::set< enum einsim::ecc_scheme >
            (
                { 
                      einsim::ES_REPETITION_T1
                    , einsim::ES_REPETITION_T2
                    , einsim::ES_REPETITION_T3
                    , einsim::ES_HAMMING_SEC
                    , einsim::ES_BCH_T1
                    , einsim::ES_BCH_T2, einsim::ES_BCH_T3
                }
            );
        }
        else
        {
            std::vector< std::string > scheme_list = options["ecc_scheme"].as< std::vector< std::string > >();
            for(std::string &ecc_scheme_str : scheme_list)
            {
                std::transform(ecc_scheme_str.begin(), ecc_scheme_str.end(), ecc_scheme_str.begin(), ::toupper);
                enum einsim::ecc_scheme ecc_scheme = einsim::ES_UNKNOWN;
                     if(ecc_scheme_str == "REP_T1") ecc_scheme = einsim::ES_REPETITION_T1;            
                else if(ecc_scheme_str == "REP_T2") ecc_scheme = einsim::ES_REPETITION_T2;            
                else if(ecc_scheme_str == "REP_T3") ecc_scheme = einsim::ES_REPETITION_T3;            
                else if(ecc_scheme_str == "HSC") ecc_scheme = einsim::ES_HAMMING_SEC;          
                else if(ecc_scheme_str == "BCH_T1") ecc_scheme = einsim::ES_BCH_T1;     
                else if(ecc_scheme_str == "BCH_T2") ecc_scheme = einsim::ES_BCH_T2;     
                else if(ecc_scheme_str == "BCH_T3") ecc_scheme = einsim::ES_BCH_T3;
                else 
                {
                    printf("[ERROR] Invalid ecc scheme requested: %s\n", ecc_scheme_str.c_str());
                    std::cout << option_parser.help({"", "Simulation"}) << std::endl;
                    return -1;
                }     
                ecc_schemes.insert(ecc_scheme);
            }
        }

        // observables - we currently only implement one, so we will always use it
        std::set< enum einsim::observable > measurements= {einsim::OBS_N_ERRORS_PER_BURST};

        /************************************************************************
         * launch the desired configuration
         ************************************************************************/
        uint64_t n_configs = data_patterns.size() * error_distributions.size() * rbers.size() * n_ecc_data_bits.size() * permutations.size() * ecc_schemes.size();
        printf("[INFO] testing %" PRIu64 " configurtions:\n", n_configs);
        std::cout << "[INFO]    " << 1                          << " burst_length_bits:   [ " << bl << " ]" << std::endl;
        std::cout << "[INFO]    " << data_patterns.size()       << " data_patterns:       [ "; for(const auto &i : data_patterns)       std::cout << einsim::enum_to_str_data_pattern(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << error_distributions.size() << " error_distributions: [ "; for(const auto &i : error_distributions) std::cout << einsim::enum_to_str_error_distribution(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << rbers.size()               << " rbers:               [ "; for(const auto &i : rbers)               std::cout << i << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << n_ecc_data_bits.size()     << " data-words:          [ "; for(const auto &i : n_ecc_data_bits)     std::cout << i << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << permutations.size()        << " permutations:        [ "; for(const auto &i : permutations)        std::cout << i << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << ecc_schemes.size()         << " ecc_schemes:         [ "; for(const auto &i : ecc_schemes)         std::cout << einsim::enum_to_str_ecc_scheme(i) << ' '; std::cout << "]" << std::endl;
        std::cout << "[INFO]    " << measurements.size()        << " observables:         [ "; for(const auto &i : measurements)        std::cout << einsim::enum_to_str_observable(i) << ' '; std::cout << "]" << std::endl;
        
        simulate
        (
              n_worker_threads
            , n_words_to_simulate
            , data_patterns
            , error_distributions
            , true_anti_cell_distributions
            , bl
            , rbers
            , n_ecc_data_bits
            , permutations
            , ecc_schemes
            , measurements
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
