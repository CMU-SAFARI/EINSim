/**
 * @file simulate.h
 *
 * @brief Core parameterized simulation routine for simulating ECC words
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef SIMULATE_H
#define SIMULATE_H

#include <set>

/* project includes */
#include "word_generator.h"
#include "observable.h"
#include "error_injector.h"
#include "ecc_code.h"
#include "supporting_routines.h"

/**
 * @brief simulates a single burst n_words_to_simulate times 
 * 
 * this function takes all the EINSim parameters and simulates 
 * the life of a dataword throughout the entire ECC process:
 * 
 * w_i -> Fenc -> C_i -> inject -> C_i' -> Fdec -> w_i' -> n_errs
 * 
 * this fn is the minimum unit of parallelism - one of these per pool entry
 * 
 * @param tid worker thread ID - useful for debugging and synchronization
 * @param ec pointer to the ECC code to simulate
 * @param n_words_to_simulate total number of bursts to simulate
 * @param burst_length_bits simulated burst length that is appropriately subdivided into ECC-code-sized pieces 
 * @param ed error distribution that specifies how errors will be injected into each burst
 * @param cd true-/anti-cell distribution to simulate
 * @param dp data pattern to simulate
 * @param rber raw bit error rate (RBER) (i.e., pre-correction error rate) to simulate
 * @param measurements list of observables to measure and output from each simulation result
 */
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
);

/**
 * @brief primary ECC simulation loop that sweeps over the provided parameters 
 * 
 * The primary simulation loop that simulates all combinations of the input
 * parameters and outputs the distribution of simulated measurements. Iterates
 * over all combinations of the provided parameters.
 * 
 * @param n_worker_threads number of worker threads to parallelize over
 * @param n_words_to_simulate  number of words to simulate for each parameter combination
 * @param data_patterns list of data patterns to test
 * @param error_distributions list of error distributions to test                           
 * @param ta_cell_distributions list of true-/anti-cell distributions to test                             
 * @param burst_length_bits list of burst lengths to test                         
 * @param rbers list of raw bit error rates (RBER) (i.e., pre-correction error rates) to test             
 * @param n_ecc_data_bits list of dataword lengths to test                       
 * @param permutations list of permutations to test                    
 * @param ecc_schemes list of enumerated ECC schemes to test                   
 * @param measurements list of observable measurements to test                    
 */
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
);

#endif /* SIMULATE_H */