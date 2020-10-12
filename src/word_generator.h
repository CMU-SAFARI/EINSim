/**
 * @file word_generator.h
 *
 * @brief Module for generating words according to specific data patterns
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef WORD_GENERATOR_H
#define WORD_GENERATOR_H

/* libraries */
#include "Eigen/Eigen"

/* project includes */
#include "supporting_routines.h"

namespace einsim
{
    /** @brief enumerates different data patterns */
    enum data_pattern
    {
          DP_RANDOM /**< each cell is randomly set to 1 or 0 */
        , DP_ALL_ONES /**< 0xFF */
        , DP_CHARGED /**< programs all cells to their charged state - every cell can experience data-retention error */
        , DP_CUSTOM /**< uses a fixed, custom pattern  */

        , DP_UNKNOWN
    };

    /** @brief enumerates different true/anti-cell distribution that we can simulate for error injection */
    enum true_anti_cell_distribution
    {
          CD_ALL_TRUE_OR_ALL_ANTI /**< 50% true-cells and 50% anti-cells, and the each burst is entirely true- or anti-cell */
        , CD_ALL_TRUE /**< 100% true cells */
        , CD_ALL_ANTI /**< 100% anti cells */
        , CD_COLSTRIPE_T /**< alternating cells starting with TRUE */
        , CD_COLSTRIPE_A /**< alternating cells starting with ANTI */

        , CD_UNKNOWN
    };

    /** @brief state representation for the true-/anti-cell layout of a given word */
    enum true_anti_cell_state
    {
          TACS_ALL_TRUE /**< the word is 100% true cells */
        , TACS_ALL_ANTI /**< the word is 100% anti cells */
        , TACS_ALT_T /**< alternating true-/anti-cells per bit starting with T */
        , TACS_ALT_A /**< alternating true-/anti-cells per bit starting with F  */

        , TACS_UNKNOWN
    };

    /** @brief state representation for the word-to-burst mapping */
    enum word_to_burst_mapping
    {
          W2BM_BLOCKS /**< the burst consists of blocks of ECC words */
        
        , W2BM_UNKNOWN
    };

    /* routines for converting between the string and enum representations of data patterns */
    enum einsim::data_pattern str_to_enum_data_pattern(const std::string &str); /**< @brief converts data_pattern string to enumeration */
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > custom_dp_to_vector(const std::string &str); /**< @brief converts a custom data pattern into an Eigen vector */
    std::string custom_dp_to_str(Eigen::Matrix< ET, Eigen::Dynamic, 1 > dp); /**< @brief converts a custom data pattern into an String representation suitable for output */
    std::string enum_to_str_data_pattern(enum einsim::data_pattern dp); /**< @brief converts data_pattern enumeration to string */
    std::string get_all_possible_data_patterns(void); /**< @brief returns a comma-separated string of all data_pattern values */
    
    /* routines for converting between the string and enum representations of true-/anti-cell distributions */
    enum einsim::true_anti_cell_distribution str_to_enum_true_anti_cell_distribution(const std::string &str); /**< @brief converts true_anti_cell_distribution string to enumeration */
    std::string enum_to_str_true_anti_cell_distribution(enum einsim::true_anti_cell_distribution cd); /**< @brief converts true_anti_cell_distribution enumeration to string */
    std::string get_all_possible_true_anti_cell_distributions(void); /**< @brief returns a comma-separated string of all true_anti_cell_distribution values */
    
    /* routines for converting between the string and enum representations of true-/anti-cell distribution state */
    enum einsim::true_anti_cell_state str_to_enum_true_anti_cell_state(const std::string &str); /**< @brief converts true_anti_cell_state string to enumeration */
    std::string enum_to_str_true_anti_cell_state(enum einsim::true_anti_cell_state cd); /**< @brief converts true_anti_cell_state enumeration to string */
    std::string get_all_possible_true_anti_cell_states(void); /**< @brief returns a comma-separated string of all true_anti_cell_state values */

    /* routines for converting between the string and enum representations of word to burst mapping */
    enum einsim::word_to_burst_mapping str_to_enum_word_to_burst_mapping(const std::string &str); /**< @brief converts word_to_burst_mapping string to enumeration */
    std::string enum_to_str_word_to_burst_mapping(enum einsim::word_to_burst_mapping w2bm); /**< @brief converts word_to_burst_mapping enumeration to string */
    std::string get_all_possible_word_to_burst_mappings(void); /**< @brief returns a comma-separated string of all word_to_burst_mapping values */

    /**
     * @brief initializes a word with the requested data pattern
     * 
     * @param word word to initialize with the requested data pattern
     * @param dp data pattern to use
     * @param custom_dp custom data pattern in the event that DP_CUSTOM is used
     * @param cd the distribution of true-/anti-cells throughout the device to sample from
     * 
     * @return returns the true-/anti-cell distribution of the genreated word
     */
    einsim::true_anti_cell_state generate_word(
          Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
        , const enum einsim::data_pattern dp
        , const Eigen::Matrix< ET, Eigen::Dynamic, 1 > &custom_dp
        , const enum einsim::true_anti_cell_distribution cd);
}

#endif /* WORD_GENERATOR_H */