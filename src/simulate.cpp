/**
 * @file simulate.cpp
 *
 * @brief Definitions for core simulation mode functions
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
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>
#include <stdarg.h>

/* libraries */
#include "libtp.h"
#include "Eigen/Eigen"

/* project includes */
#include "word_generator.h"
#include "observable.h"
#include "error_model.h"
#include "ecc_code.h"
#include "supporting_routines.h"
#include "codes/repetition_code.h"
#include "codes/bch_code.h"
#include "codes/hamming_code.h"

void simulate_burst
(
      const int tid
    , const einsim::ecc_code *ec // ECC code to simulate
    , const uint64_t n_bursts_to_simulate
    , const int burst_length_bits
    , const enum einsim::word_to_burst_mapping w2b_map
    , const std::vector< einsim::error_model_descriptor * > &emd
    , const enum einsim::true_anti_cell_distribution cd
    , const enum einsim::data_pattern dp
    , const Eigen::Matrix< ET, Eigen::Dynamic, 1 > &custom_dp
    , const std::set< enum einsim::observable > &observables
)
{
    int pad_size = (ec->get_n_data_bits() - (burst_length_bits % ec->get_n_data_bits())) % ec->get_n_data_bits();
    int n_ecc_words_per_burst = (burst_length_bits / ec->get_n_data_bits()) + !!pad_size;
    int burst_codeword_length_bits = ec->get_n_code_bits() * n_ecc_words_per_burst; // NOTE: we should subtract pad_size here, but our ECC
                                                                                    // model requires passing all the data bits. Instead,
                                                                                    // we pass the pad bits to the outfile so it can be analyzed as necessary

    // declare intermediate words
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > burst_word;
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > dataword[n_ecc_words_per_burst];
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > codeword[n_ecc_words_per_burst];
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > burst_codeword;
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > burst_codeword_p;
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > codeword_p[n_ecc_words_per_burst];
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > dataword_p[n_ecc_words_per_burst];
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > burst_word_p;

    // size intermediate words
    burst_word = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(burst_length_bits);
    for(int i = 0; i < n_ecc_words_per_burst; i++)
    {
        dataword[i] = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(ec->get_n_data_bits());
        codeword[i] = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(ec->get_n_code_bits());
    }
    burst_codeword = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(burst_codeword_length_bits);
    burst_codeword_p = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(burst_codeword_length_bits);
    for(int i = 0; i < n_ecc_words_per_burst; i++)
    {
        codeword_p[i] = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(ec->get_n_code_bits());
        dataword_p[i] = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(ec->get_n_data_bits());
    }
    burst_word_p = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(burst_length_bits);

    // begin simulation
    std::map< uint32_t /* nerrs */, std::pair< uint64_t /* rber */, uint64_t /* uber */ > > accumulator;
    std::vector< uint64_t > per_bit_histogram_burst_data(burst_length_bits);
    std::vector< uint64_t > per_bit_histogram_burst_code(burst_codeword_length_bits);
    for(uint64_t n_bursts_simulated = 0; n_bursts_simulated < n_bursts_to_simulate; n_bursts_simulated++)
    {       
        // Generate a burst
        einsim::true_anti_cell_state tacs = einsim::generate_word(burst_word, dp, custom_dp, cd);

        // Split the burst into datawords (try for equal splitting if at all possible)
        for(int i = 0; i < n_ecc_words_per_burst; i++)
        {
            // extract the dataword (padding if necessary)
            if(pad_size && i == n_ecc_words_per_burst - 1)
            {
                dataword[i].segment(0, ec->get_n_data_bits() - pad_size) = 
                    burst_word.segment(ec->get_n_data_bits() * i, ec->get_n_data_bits() - pad_size);
                for(int j = ec->get_n_data_bits() - pad_size; j < dataword[i].size(); j++)
                    dataword[i][j] = 0; // pad
            }
            else
                dataword[i] = burst_word.segment(ec->get_n_data_bits() * i, ec->get_n_data_bits());
            // std::cout << "Dataword[" << i << "]: " << dataword[i].transpose() << std::endl;
        }

        // encode each dataword
        for(int i = 0; i < n_ecc_words_per_burst; i++)
            codeword[i] = ec->encode(dataword[i]);

        // assemble the burst codeword
        for(int i = 0; i < n_ecc_words_per_burst; i++)
            burst_codeword.segment(ec->get_n_code_bits() * i, ec->get_n_code_bits()) = codeword[i];
        
        // induce errors in the burst_codeword
        burst_codeword_p = burst_codeword;
        einsim::inject(burst_codeword_p, dp, tacs, emd);

        // extract the individual codewords from the burst codeword
        for(int i = 0; i < n_ecc_words_per_burst; i++)
            codeword_p[i] = burst_codeword_p.segment(ec->get_n_code_bits() * i, ec->get_n_code_bits());

        // correct the codewords
        for(int i = 0; i < n_ecc_words_per_burst; i++)
        {
            dataword_p[i] = ec->decode(codeword_p[i]);
            
            // quick sanity check
            int nerrs_induced = hamming_distance(codeword[i], codeword_p[i]);
            int nerrs_observed = hamming_distance(dataword[i], dataword_p[i]);
            bool fully_correctable = ec->correction_capability() >= nerrs_induced;
            if(fully_correctable && nerrs_observed) // sanity check
            {
                printf_both("[ERROR]: observed %d errors when %d errors induced and %d correctable (%s)\n", nerrs_observed, nerrs_induced, ec->correction_capability(), ec->name().c_str());
                std::stringstream ss;
                ss << "code_sent: " << codeword[i].transpose() << ", code_rcvd: " << codeword_p[i].transpose()
                    << ", data_sent: " << dataword[i].transpose() << ", data_rcvd: " << dataword_p[i].transpose();
                printf_both("[ERROR] %s\n", ss.str().c_str());
                assert(0);
            }
        }

        // reassemble the burst
        for(int i = 0; i < n_ecc_words_per_burst; i++)
        {
            // handle necessary padding on the last word
            if(pad_size && i == n_ecc_words_per_burst - 1)
            {
                burst_word_p.segment(ec->get_n_data_bits() * i, ec->get_n_data_bits() - pad_size) = 
                    dataword_p[i].segment(0, ec->get_n_data_bits() - pad_size);
            }
            else
                burst_word_p.segment(ec->get_n_data_bits() * i, ec->get_n_data_bits()) = dataword_p[i];
        }

        // print statii
        if(g_verbosity >= 2)
        {
            std::cout << "Burst    : " << burst_word.transpose() << std::endl;
            std::cout << "Burst cw : " << burst_codeword.transpose() << std::endl;
            std::cout << "Burst cw': " << burst_codeword_p.transpose() << std::endl;
            std::cout << "Burst'   : " << burst_word_p.transpose() << std::endl;
            for(int i = 0; i < n_ecc_words_per_burst; i++)
            {
                std::cout << "    w[" << i << "] : " << dataword[i].transpose() << std::endl;
                std::cout << "    c[" << i << "] : " << codeword[i].transpose() << std::endl;
                std::cout << "    c[" << i << "]': " << codeword_p[i].transpose() << std::endl;
                std::cout << "    w[" << i << "]': " << dataword_p[i].transpose() << std::endl;
            }
        }

        // account for the RBER
        int nerrs_in_burst_codeword = hamming_distance(burst_codeword, burst_codeword_p);
        if(accumulator.count(nerrs_in_burst_codeword) == 0)
            accumulator[nerrs_in_burst_codeword] = std::make_pair(0, 0);
        std::get< 0 >(accumulator[nerrs_in_burst_codeword]) += 1;
        
        // account for the UBER
        int nerrs_in_burst = hamming_distance(burst_word, burst_word_p);
        if(accumulator.count(nerrs_in_burst) == 0)
            accumulator[nerrs_in_burst] = std::make_pair(0, 0);
        std::get< 1 >(accumulator[nerrs_in_burst]) += 1;

        // fill in the per-bit histograms
        for(int i = 0; i < burst_length_bits; i++)
            if(burst_word[i] ^ burst_word_p[i])
                per_bit_histogram_burst_data[i]++;
        for(int i = 0; i < burst_codeword_length_bits; i++)
            if(burst_codeword[i] ^ burst_codeword_p[i])
                per_bit_histogram_burst_code[i]++;
    }

    // as of now, we only have one observable
    for(const auto &obs : observables)
    {
        // prepare a line of output representing the results of this simulation
        std::stringstream ss;
        ss << "[DATA] uid:" << ec->get_uid()
            << " nw:" << n_bursts_to_simulate
            << " bl:" << burst_length_bits
            << " bcl:" << burst_codeword_length_bits
            << " ps:" << pad_size
            << " em:" << einsim::error_model_descriptor_vec_to_str(emd)
            << " cd:" << einsim::enum_to_str_true_anti_cell_distribution(cd)
            << " dp:" << einsim::enum_to_str_data_pattern(dp);
            if(dp == einsim::DP_CUSTOM)
                ss << " cdp:" << einsim::custom_dp_to_str(custom_dp);
            ss << " obs:" << einsim::enum_to_str_observable(obs);
        switch(obs)
        {
            case einsim::OBS_N_ERRORS_PER_BURST:
            {
                ss << " [ ";
                for(const auto &accum_it : accumulator)
                {
                    uint32_t nerrs = accum_it.first;
                    uint64_t re_count = std::get< 0 >(accum_it.second);
                    uint64_t ue_count = std::get< 1 >(accum_it.second);
                    ss << nerrs << ':' << re_count << ':' << ue_count << ' ';
                }
                ss << "]";
                break;
            }

            case einsim::OBS_PER_BIT_ERROR_COUNT:
            {
                ss << " [ ";
                for(int i = 0; i < burst_length_bits; i++)
                    ss << per_bit_histogram_burst_data[i] << ' ';
                ss << ": ";
                for(int i = 0; i < burst_codeword_length_bits; i++)
                    ss << per_bit_histogram_burst_code[i] << ' ';
                ss << "]";
                break;
            }

            default:
                printf_both("[ERROR] unsupported observable: %d\n", (int)obs);
                assert(0);
        }

        // print to the output file (either stdout or the file)
        // printf_both("%s\n", ss.str().c_str());
        fprintf(g_output_file, "%s\n", ss.str().c_str());
        if(g_verbosity > 1)
            printf("%s\n", ss.str().c_str());
    } 

}

void simulate
(
      const int n_worker_threads
    , const uint64_t n_bursts_to_simulate
    , const uint64_t n_bursts_per_job
    , const std::set< int > &burst_lengths_nbits
    , const std::set< enum einsim::word_to_burst_mapping > &w2b_mappings
    , const std::vector< enum einsim::data_pattern > &data_patterns
    , const std::vector< Eigen::Matrix< ET, Eigen::Dynamic, 1 > > &custom_patterns
    , const std::vector< std::vector< einsim::error_model_descriptor * > > &error_models
    , const std::set< enum einsim::true_anti_cell_distribution > &ta_cell_distributions
    , const std::set< enum einsim::observable > &observables
    , const std::vector< einsim::ecc_code * > &ecc_schemes
)
{
    thread_pool tp(n_worker_threads);
    tp.start();

    // std::set< int > data_bits = {16, 32, 64, 128, 256, 512, 1024};
    // std::map< int /* n_data_bits */, std::vector< einsim::ecc_code * > > ecc_types;

    // print out all ECC schemes so that we can map UIDs to them
    for(const einsim::ecc_code *ec : ecc_schemes)
    {
        std::string json;
        if(ec->to_json(json))
        {
            std::cout << "[ERROR] could not convert ECC scheme " << ec->name() << " to JSON" << std::endl;
            assert(0 && "failed to convert ECC scheme to JSON format");
        }
        json.erase(std::remove_if(json.begin(), json.end(), ::isspace), json.end());
        fprintf(g_output_file, "[ECC] %s\n", json.c_str());
    }

    // print out what is happening
    if(g_verbosity > 0)
        printf("[INFO] Starting ECC simulations\n");
    fprintf(g_output_file, "[INFO] Starting ECC simulations\n");

    // sweep the requested parameters and begin simulations
    int custom_dp_idx = 0;
    for(const enum einsim::data_pattern dp : data_patterns)
    {
        for(const int bl : burst_lengths_nbits)
        {
            // determine the custom data pattern of length bl
            Eigen::Matrix< ET, Eigen::Dynamic, 1 > custom_dp = Eigen::Matrix< ET, Eigen::Dynamic, 1 >::Zero(bl);
            if(dp == einsim::DP_CUSTOM)
            {
                const Eigen::Matrix< ET, Eigen::Dynamic, 1 > &custom_dp_raw = custom_patterns.at(custom_dp_idx);
                if(custom_dp_raw.rows() != bl)
                {
                    printf_both("[ERROR] providing a custom DP must exactly match the burst length. bl: %d, dp_length: %d\n", bl, custom_dp_raw.rows());
                    assert(0 && "custom data pattern must be exactly the length of a burst");  
                } 
                // custom_dp.block(custom_dp.rows() - custom_dp_raw.rows(), 0, custom_dp_raw.rows(), 1) = custom_dp_raw;
                custom_dp = custom_dp_raw;
                assert(custom_dp.rows() == bl && custom_dp.cols() == 1);
            }

            for(const enum einsim::word_to_burst_mapping &w2b_map : w2b_mappings)
            {
                for(const std::vector< einsim::error_model_descriptor *> em : error_models)
                {
                    for(const enum einsim::true_anti_cell_distribution cd : ta_cell_distributions)
                    {
                        for(const einsim::ecc_code *ec : ecc_schemes)
                        {
                            if(!(em.size() == 1 || (int)em.size() == ec->get_n_code_bits()))
                            {
                                std::cout << "[ERROR] each error model size must be 1 or match the codeword size. em.size(): " << em.size() << " codeword size: " << ec->get_n_code_bits() << std::endl;
                                assert(0 && "error model size must be 1 or match the codeword size");  
                            }

                            uint64_t n_bursts_remaining = n_bursts_to_simulate;
                            while(n_bursts_remaining)
                            {
                                uint64_t n_bursts_this_job = std::min(n_bursts_remaining, n_bursts_per_job);
                                tp.add(simulate_burst, 0 /* priority */, ec, n_bursts_this_job
                                    , bl, w2b_map, em, cd, dp, custom_dp, observables); /* bursts are accumulated */
                                n_bursts_remaining -= n_bursts_this_job;
                                
                                // avoid overflowing memory with the number of requested jobs
                                while(tp.get_n_jobs_outstanding() > 1000000) // 1GiB worth
                                    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                            }
                        }
                    }
                }
            }
        }

        // the custom DP array is indexed 0..N-1 only for custom DPs
        if(dp == einsim::DP_CUSTOM)
            custom_dp_idx++;
    }

    // wait for outstanding jobs to complete
    auto start_time = std::chrono::system_clock::now();
    int n_outstanding;
    while((n_outstanding = tp.get_n_jobs_outstanding()))
    {   
        auto elapsed_time = std::chrono::system_clock::now() - start_time;
        int n_total_jobs = n_outstanding + tp.get_n_jobs_completed(); // a bit of harmless TOCTOU
        float ms_per_job = std::chrono::duration_cast< std::chrono::milliseconds >(elapsed_time).count()  / (float)tp.get_n_jobs_completed();
        float eta_ms = ms_per_job * n_outstanding;
        printf("Jobs remaining: %d/%d (ETA: %d:%02d:%02d)\n", n_outstanding, n_total_jobs
            , (int)(eta_ms / (1000 * 60 * 60))
            , (int)(eta_ms / (1000 * 60)) % 60
            , (int)(eta_ms / (1000)) % 60);
        std::this_thread::sleep_for(std::chrono::milliseconds(n_outstanding <= 8 ? 64 : 2000));
    }
    tp.wait();
    tp.reset_stats();
}
