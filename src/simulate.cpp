/**
 * @file simulate.cpp
 *
 * @brief Definitions for core simulation mode functions
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <cstdio>
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

/* project includes */
#include "word_generator.h"
#include "observable.h"
#include "error_injector.h"
#include "ecc_code.h"
#include "supporting_routines.h"
#include "codes/repetition_code.h"
#include "codes/bch_code.h"
#include "codes/hamming_code.h"

void simulate_burst
(
	  int tid
	, einsim::ecc_code *ec // ECC code to simulate
	, int n_words_to_simulate
	, int burst_length_bits
	, enum einsim::error_distribution ed
	, enum einsim::true_anti_cell_distribution cd
	, enum einsim::data_pattern dp
	, float rber
    , std::set< enum einsim::observable > measurements
)
{
    std::map< uint32_t /* nerrs */, std::pair< uint64_t /* rber */, uint64_t /* uber */ > > accumulator;
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
    for(int n_words_simulated = 0; n_words_simulated < n_words_to_simulate; n_words_simulated++)
    {       
        // Generate a burst
        einsim::generate_word(burst_word, dp);

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
        einsim::inject(burst_codeword_p, ed, cd, dp, rber);

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
	            exit(-1);
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
    }

    // as of now, we only have one observable
    assert(measurements.size() == 1 && *measurements.begin() == einsim::OBS_N_ERRORS_PER_BURST);

    // prepare a line of output representing the results of this simulation
    std::stringstream ss;
    ss << "[DATA] " << ec->name_short() // also prints permutation
        << " rber:" << rber
        << " bl:" << burst_length_bits
        << " bcl:" << burst_codeword_length_bits
        << " ps:" << pad_size
        << " ed:" << einsim::enum_to_str_error_distribution(ed)
        << " cd:" << einsim::enum_to_str_true_anti_cell_distribution(cd)
        << " dp:" << einsim::enum_to_str_data_pattern(dp) 
        << " [ ";
    for(const auto &accum_it : accumulator)
    {
        uint32_t nerrs = accum_it.first;
        uint64_t re_count = std::get< 0 >(accum_it.second);
        uint64_t ue_count = std::get< 1 >(accum_it.second);
        ss << nerrs << ':' << re_count << ':' << ue_count << ' ';
    }
    ss << "]";

    // print to the output file (either stdout or the file)
    // printf_both("%s\n", ss.str().c_str());
    fprintf(g_output_file, "%s\n", ss.str().c_str());
    if(g_verbosity > 1)
        printf("%s\n", ss.str().c_str());
}

void simulate
(
      int n_worker_threads
    , int n_words_to_simulate
    , std::set< enum einsim::data_pattern > data_patterns
    , std::set< enum einsim::error_distribution > error_distributions
    , std::set< enum einsim::true_anti_cell_distribution > ta_cell_distributions
    , int burst_length_bits
    , std::set< float > rbers
    , std::set< int > n_ecc_data_bits
    , std::set< int > permutations
    , std::set< enum einsim::ecc_scheme > ecc_schemes
    , std::set< enum einsim::observable > measurements
)
{
    thread_pool tp(n_worker_threads);
    tp.start();

    // std::set< int > data_bits = {16, 32, 64, 128, 256, 512, 1024};
    for(int ecc_permutation : permutations)
    {
        std::map< int /* n_data_bits */, std::vector< einsim::ecc_code * > > ecc_types;

        /* prepare the current ECC types */
        if(g_verbosity > 0)
            printf("Preparing ECC codes...\n");
        fprintf(g_output_file, "Preparing ECC codes...\n");

        // instantiate the different ECC types for the data bit sizes and ECC schemes requested
        for(int n_data_bits : n_ecc_data_bits)
        {
            for(enum einsim::ecc_scheme scheme : ecc_schemes)
            {
                switch(scheme)
                {
                    case einsim::ES_REPETITION_T1:
                    case einsim::ES_REPETITION_T2:
                    case einsim::ES_REPETITION_T3:
                    {
                        int correction_capability = (scheme == einsim::ES_REPETITION_T1) ? 3 : ((scheme == einsim::ES_REPETITION_T2) ? 5 : 7);
                        ecc_types[n_data_bits].push_back(new einsim::repetition(ecc_permutation, n_data_bits, correction_capability));
                        break;
                    }
                    case einsim::ES_HAMMING_SEC:
                    {
                        ecc_types[n_data_bits].push_back(new einsim::hamming(ecc_permutation, n_data_bits));
                        break;
                    }
                    case einsim::ES_BCH_T1:
                    case einsim::ES_BCH_T2:
                    case einsim::ES_BCH_T3:
                    {
                        int correction_capability = (scheme == einsim::ES_BCH_T1) ? 1 : ((scheme == einsim::ES_BCH_T2) ? 2 : 3);
                        einsim::bch *bch_code = new einsim::bch(ecc_permutation, n_data_bits, correction_capability);
                        if(bch_code->ready())
                            ecc_types[n_data_bits].push_back(bch_code);
                        else
                        {
                            printf_both("[ERROR] no such code exists for the requested parameters!\n");
                            exit(-1);
                        }
                        break;
                    }
                    default:
                        printf("[ERROR] unknown/invalid ECC scheme: %d\n", scheme);
                        exit(-1);
                }
            }
        }
        
        if(g_verbosity > 0)
            printf("Starting ECC simulations\n");
        fprintf(g_output_file, "Starting ECC simulations\n");

        // sweep the requested parameters and begin simulations
        for(enum einsim::data_pattern dp : data_patterns)
        {
            for(enum einsim::error_distribution ed : error_distributions)
            {
                for(enum einsim::true_anti_cell_distribution cd : ta_cell_distributions)
                {
                    for(int n_data_bits : n_ecc_data_bits)
                    {
                        for(unsigned et_idx = 0; et_idx < ecc_types[n_data_bits].size(); et_idx++)
                        {
                            einsim::ecc_code *ec = ecc_types[n_data_bits][et_idx];
                            
                            for(float rber : rbers)
                            {
                                // return_values[dp][em][n_data_bits][et_idx] = accum;

                                // some simple load balancing - some ECCs are more expensive than others :)
                                int n_jobs = 1;
                                int n_words_per_job = n_words_to_simulate;
                                int max_n_words_per_job = 10000; // sort of a parallelization factor
                                if(n_words_to_simulate > max_n_words_per_job)
                                {
                                    n_jobs = 1 + (n_words_to_simulate / max_n_words_per_job); // MIGHT overdo the total number!!
                                    n_words_per_job = max_n_words_per_job;
                                }
                                for(int n = 0; n < n_jobs; n++)
                                    tp.add(simulate_burst, 0 /* priority */, ec, n_words_per_job, burst_length_bits, ed, cd, dp, rber, measurements); /* bursts are accumulated */
                            }
                        }
                    }
                }
            }
        }

        // wait for outstanding jobs to complete
        while(tp.get_n_jobs_outstanding())
        {   
            int n_total_jobs = tp.get_n_jobs_outstanding() + tp.get_n_jobs_completed();
            printf("Jobs remaining: %d/%d\n", tp.get_n_jobs_outstanding(), n_total_jobs);
            std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        }
        tp.wait();
        tp.reset_stats();
    }
}