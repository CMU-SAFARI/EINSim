/**
 * @file error_injector.h
 *
 * @brief Data structures and routines used for injecting parameterizable errors
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef ERROR_INJECTOR_H
#define ERROR_INJECTOR_H

#include <string>
#include <algorithm>
#include <random>

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

		, DP_UNKNOWN
	};

	/** @brief enumerates different error distributions that the error injector can use */
	enum error_distribution
	{
		  ED_UNIFORM_RANDOM /**< uniform-random error injection */

		, ED_UNKNOWN
	};

	/** @brief enumerates different true/anti-cell distribution that we can simulate for error injection */
	enum true_anti_cell_distribution
	{
		  CD_ALL_TRUE_OR_ALL_ANTI /**< 50% true-cells and 50% anti-cells, and the each burst is entirely true- or anti-cell */
		, CD_ALL_TRUE /**< 100% true cells */
		, CD_ALL_ANTI /**< 100% anti cells */

		, CD_UNKNOWN
	};

	/* routines for converting between the string and enum representations of data patterns */
	enum einsim::data_pattern str_to_enum_data_pattern(std::string str); /**< @brief converts data_pattern string to enumeration */
	std::string enum_to_str_data_pattern(enum einsim::data_pattern dp); /**< @brief converts data_pattern enumeration to string */
	std::string get_all_possible_data_patterns(void); /**< @brief returns a comma-separated string of all data_pattern values */
	
	/* routines for converting between the string and enum representations of error distributions */
	enum einsim::error_distribution str_to_enum_error_distribution(std::string str); /**< @brief converts error_distribution string to enumeration */
	std::string enum_to_str_error_distribution(enum einsim::error_distribution ed); /**< @brief converts error_distribution enumeration to string */
	std::string get_all_possible_error_distributions(void); /**< @brief returns a comma-separated string of all error_distribution values */
	
	/* routines for converting between the string and enum representations of true-/anti-cell distributions */
	enum einsim::true_anti_cell_distribution str_to_enum_true_anti_cell_distribution(std::string str); /**< @brief converts true_anti_cell_distribution string to enumeration */
	std::string enum_to_str_true_anti_cell_distribution(enum einsim::true_anti_cell_distribution cd); /**< @brief converts true_anti_cell_distribution enumeration to string */
	std::string get_all_possible_true_anti_cell_distributions(void); /**< @brief returns a comma-separated string of all true_anti_cell_distribution values */

	/**
	 * @brief probabilistically injects errors into a word with the given parameters
	 * 
	 * Injects errors according to the specified error distribution. This is
	 * the primary error injection routine that simulates pre-correction error
	 * characteristics.
	 * 
	 * @param word the word in which to inject errors
	 * @param ed the error distribution 
	 * @param cd the true-/anti-cell distribution to assume when injecting errors
	 * @param dp the programmed data pattern to assume when injecting errors
	 * @param rber the raw bit error rate (RBER) (i.e., pre-correction error rate) to target for error injection
	 */
	void inject(
		  Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
		, enum error_distribution ed
		, enum true_anti_cell_distribution cd
		, enum data_pattern dp
		, float rber
	);
	
	/**
	 * @brief injects a precise number of errors into a word with the given parameters
	 * 
	 * Generally used for testing and debugging, this function allows
	 * injecting a precise number of errors into a word, though the actual
	 * location of injected errors will vary according to the provided
	 * parameters.
	 * 
	 * @param word the word in which to inject errors
	 * @param ed the error distribution 
	 * @param cd the true-/anti-cell distribution to assume when injecting errors
	 * @param dp the programmed data pattern to assume when injecting errors
	 * @param n_errors the number of errors to inject
	 */
	void inject_n(
		  Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
		, enum error_distribution ed
		, enum true_anti_cell_distribution cd
		, enum data_pattern dp
		, int n_errors
	);
}

#endif /* ERROR_INJECTOR_H */