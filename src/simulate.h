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
#include "error_model.h"
#include "ecc_code.h"
#include "supporting_routines.h"

/**
 * @brief simulates a single burst n_bursts_to_simulate times 
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
 * @param n_bursts_to_simulate total number of bursts to simulate
 * @param burst_length_bits simulated burst length that is appropriately subdivided into ECC-code-sized pieces 
 * @param w2b_map mapping by which data is mapped to a burst 
 * @param emd error model that specifies how errors will be injected into each burst
 * @param cd true-/anti-cell distribution to simulate
 * @param dp data pattern to simulate
 * @param custom_dp custom data pattern to use in the event that it's been specified
 * @param observables list of observables to measure and output from each simulation result
 */
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
);

/**
 * @brief primary ECC simulation loop that sweeps over the provided parameters
 * 
 * The primary simulation loop that simulates all combinations of the input
 * parameters and outputs the distribution of simulated measurements. Iterates
 * over all combinations of the provided parameters.
 * 
 * @param n_worker_threads number of worker threads to parallelize over
 * @param n_bursts_to_simulate number of words to simulate for each parameter combination
 * @param n_bursts_per_job number of bursts to simulate per group
 * @param burst_lengths_nbits list of burst lengths to test
 * @param w2b_mappings list of mappings by which data is mapped to a burst 
 * @param data_patterns list of data patterns to test
 * @param custom_patterns list of custom data patterns
 * @param error_models list of error model descriptors to test
 * @param ta_cell_distributions list of true-/anti-cell distributions to test
 * @param observables list of observable observables to measure
 * @param ecc_schemes list of ECC schemes to test
 */
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
);

#endif /* SIMULATE_H */