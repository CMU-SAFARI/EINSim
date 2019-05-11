/**
 * @file word_generator.cpp
 *
 * @brief Routines for generating words according to specific data patterns
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <string>
#include <algorithm>
#include <random>

/* libraries */
#include "Eigen/Eigen"

/* project includes */
#include "supporting_routines.h"
#include "word_generator.h"

static std::default_random_engine word_generator_rng;
static std::bernoulli_distribution bernoulli_0_5(0.5);
static bool word_generator_initialized = false;

void einsim::generate_word(Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word, enum einsim::data_pattern dp)
{
	// handle the one-time initialization of the global random number generator
	if(!word_generator_initialized)
	{
		word_generator_rng.seed(std::random_device{}());
		word_generator_initialized = true;
	}

	// write in the appropriate data pattern to the word
	switch(dp)
	{
		case einsim::DP_RANDOM:
		{
			for(int i = 0; i < word.size(); i++)
				word[i] = bernoulli_0_5(word_generator_rng) ? 0 : 1;
			break;
		}
		case einsim::DP_ALL_ONES:
		{
			word = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Ones(word.size());
			break;
		}
		case einsim::DP_CHARGED:
		{
			bool true_cells = bernoulli_0_5(word_generator_rng);
			for(int i = 0; i < word.size(); i++)
				word[i] = true_cells ? 1 : 0;
			break;
		}
		default:
		{
			printf_both("[ERROR] unknown or unsupported data pattern: %d", dp);
			exit(-1);
		}
	}
}