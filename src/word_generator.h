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
#include "error_injector.h"
#include "supporting_routines.h"

namespace einsim
{
	/**
	 * @brief initializes a word with the requested data pattern
	 * 
	 * @param word word to initialize with the requested data pattern
	 * @param dp data pattern to use
	 */
	void generate_word(Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word, enum einsim::data_pattern dp);
}

#endif /* WORD_GENERATOR_H */