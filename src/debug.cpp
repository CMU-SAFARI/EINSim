/**
 * @file debug.cpp
 *
 * @brief Provides simple simulation loops that are useful for debugging ECC implementations
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
#include "Eigen/Eigen"
#include "cxxopts.h"

/* project */
#include "word_generator.h"
#include "error_injector.h"
#include "ecc_code.h"
#include "supporting_routines.h"
#include "codes/repetition_code.h"
#include "codes/bch_code.h"
#include "codes/hamming_code.h"
#include "debug.h"

void debug_example_worker(int tid, einsim::ecc_code *ec, int n_words_to_simulate, enum einsim::data_pattern dp)
{
    std::map< uint32_t, std::map< uint32_t, uint64_t > > accumulator;

    for(int word = 0; word < n_words_to_simulate; word++)
    {       
        // Generate dataword
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_sent = 
            Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(ec->get_n_data_bits());
        einsim::generate_word(data_word_sent, dp);

        // Encode dataword into codeword
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word_sent = ec->encode(data_word_sent);

        // Corrupt codeword with N errors
        for(int nerrs_transmitted = 0; nerrs_transmitted <= ec->get_n_code_bits(); nerrs_transmitted++)
        {
            bool fully_correctable = ec->correction_capability() >= nerrs_transmitted;

            // inject precisely N errors. the DP is not strictly true, but we're just trying to test functionality here 
            Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word_recv = code_word_sent;
            inject_n(code_word_recv, einsim::ED_UNIFORM_RANDOM, einsim::CD_ALL_TRUE, einsim::DP_CHARGED, nerrs_transmitted);

            int nerrs_induced = hamming_distance(code_word_sent, code_word_recv);
            if(nerrs_induced > nerrs_transmitted)
            {
                printf_both("[ERROR] more errs induced (%d) than errs transmitted(%d) (%s)\n", nerrs_induced, nerrs_transmitted, ec->name().c_str());
                std::stringstream ss;
                ss << "code_sent: " << code_word_sent << ", code_rcvd: " << code_word_recv
                    << ", data_sent: " << data_word_sent;
                printf_both("[ERROR] %s\n", ss.str().c_str());
                exit(-1);
            }

            // Decode codeword into dataword
            Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_recv = ec->decode(code_word_recv);

            // print the messages
            // std::cout << "Dataword orig: " << data_word_sent.transpose() << std::endl;
            // std::cout << "Codeword sent: " << code_word_sent.transpose() << std::endl;
            // std::cout << "Codeword recv: " << code_word_recv.transpose() << std::endl;
            // std::cout << "Dataword recv: " << data_word_recv.transpose() << std::endl;

            int nerrs_observed = hamming_distance(data_word_sent, data_word_recv);
            if(fully_correctable && nerrs_observed) // sanity check
            {
                printf_both("[ERROR]: observed %d errors when %d induced and %d correctable (%s)\n", nerrs_observed, nerrs_transmitted, ec->correction_capability(), ec->name().c_str());
                std::stringstream ss;
                ss << "code_sent: " << code_word_sent << ", code_rcvd: " << code_word_recv
                    << ", data_sent: " << data_word_sent << ", data_rcvd: " << data_word_recv;
                printf_both("[ERROR] %s\n", ss.str().c_str());
                exit(-1);
            }
            
            if(accumulator[nerrs_transmitted].count(nerrs_observed) == 0)
                accumulator[nerrs_transmitted][nerrs_observed] = 0;
            accumulator[nerrs_transmitted][nerrs_observed]++;
       }
   }

    std::stringstream ss;
    ss << ec->name_short() << " dp:" << enum_to_str_data_pattern(dp) << " [ ";
    for(const auto &acc_it : accumulator)
    {
        uint32_t nerrs_transmitted = acc_it.first;
        for(const auto &transmit_it : acc_it.second)
        {
            uint32_t nerrs_observed = transmit_it.first;
            uint64_t count = transmit_it.second;
            ss << nerrs_transmitted << ":" << nerrs_observed << ':' << count << ' ';
        }
    }
    ss << "]";

    fprintf(g_output_file, "%s\n", ss.str().c_str());
    if(g_verbosity > 1)
        printf("%s\n", ss.str().c_str());
    return;
}

void debug_example(int n_worker_threads, int n_words_to_simulate)
{
    thread_pool tp(n_worker_threads);
    tp.start();

    std::set< int > data_bits_base = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
    std::set< int > data_bits = {};
    for(int ecc_permutation = 0;; ecc_permutation++) // will continue infinitely
    {
        std::map< int /* n_data_bits */, std::vector< einsim::ecc_code * > > ecc_types;

        // prepare some ECC codes to test
        if(g_verbosity > 0)
            printf("Preparing ECC codes...\n");
        fprintf(g_output_file, "Preparing ECC codes...\n");

        for(int n_data_bits_base : data_bits_base)
        {
            for(int i = -1; i <= 1; i++)
            {
                int n_data_bits = n_data_bits_base + i;
                if(n_data_bits < 1)
                    continue;

                ecc_types[n_data_bits].push_back(new einsim::repetition(ecc_permutation, n_data_bits, 3));
                ecc_types[n_data_bits].push_back(new einsim::hamming(ecc_permutation, n_data_bits));
                ecc_types[n_data_bits].push_back(new einsim::bch(ecc_permutation, n_data_bits, 3));
                ecc_types[n_data_bits].push_back(new einsim::bch(ecc_permutation, n_data_bits, 5));
                ecc_types[n_data_bits].push_back(new einsim::bch(ecc_permutation, n_data_bits, 7));
                ecc_types[n_data_bits].push_back(new einsim::bch(ecc_permutation, n_data_bits, 9));
                data_bits.insert(n_data_bits);
            }
        }
        
        // begin parallel simulations over the desired parameter space 
        if(g_verbosity > 0)
            printf("Starting ECC simulations\n");
        fprintf(g_output_file, "Starting ECC simulations\n");

        enum einsim::data_pattern all_dps[] = {einsim::DP_RANDOM, einsim::DP_CHARGED, einsim::DP_ALL_ONES};
        for(int dp_idx = 0; dp_idx < (int)(sizeof(all_dps) / sizeof(all_dps[0])); dp_idx++)
        {
            enum einsim::data_pattern dp = all_dps[dp_idx];
            for(int n_data_bits : data_bits)
            {
              // simulate several repetitions - the samples for identical experiments can be accumulated
                int n_repeats = 10;
                for(unsigned et_idx = 0; et_idx < ecc_types[n_data_bits].size(); et_idx++)
                {
                    einsim::ecc_code *ec = ecc_types[n_data_bits][et_idx];
                    for(int n = 0; n < n_repeats; n++)
                        tp.add(debug_example_worker, 0 /* priority */, ec, n_words_to_simulate, dp);
                }
            }
        }

        // wait for outstanding jobs to complete
        while(tp.get_n_jobs_outstanding())
        {   
            int n_total_jobs = tp.get_n_jobs_outstanding() + tp.get_n_jobs_completed();
            if(g_verbosity > 0)
            {
                printf("Testing: [%d/%d] jobs remaining     \r", tp.get_n_jobs_outstanding(), n_total_jobs);
                fflush(stdout);
            }
            fprintf(g_output_file, "Jobs remaining: %d/%d\n", tp.get_n_jobs_outstanding(), n_total_jobs);
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
        tp.wait();
        tp.reset_stats();
    }
    
}
