/**
 * @file error_injector.cpp
 *
 * @brief Definitions of routines for injecting errors into ECC words
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <iostream>
#include <string>
#include <algorithm>
#include <random>

#include "Eigen/Eigen"
#include "supporting_routines.h"
#include "error_injector.h"

/* random distribution parameters used by the error injector */
static std::default_random_engine rng;
static std::bernoulli_distribution ei_bernoulli_0_5(0.5);
static bool error_injector_initialized = false;

static std::random_device rd_inject_errs;
static std::mt19937 rng_inject_errs(rd_inject_errs());

enum einsim::data_pattern einsim::str_to_enum_data_pattern(std::string str)
{
	static std::map< std::string, enum data_pattern > str_to_enum_data_pattern_map = 
	{
		  {"RANDOM", DP_RANDOM}
		, {"ALL_ONES", DP_ALL_ONES}
		, {"CHARGED", DP_CHARGED}
	};
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
	auto search = str_to_enum_data_pattern_map.find(str);
	if(search == str_to_enum_data_pattern_map.end())
		return DP_UNKNOWN;
	else
		return search->second;
}

std::string einsim::enum_to_str_data_pattern(enum einsim::data_pattern dp)
{
	static std::map< enum data_pattern, std::string > enum_to_str_data_pattern_map = 
	{
		  {DP_RANDOM, "RANDOM"}
		, {DP_ALL_ONES, "ALL_ONES"}
		, {DP_CHARGED, "CHARGED"}
	};
	auto search = enum_to_str_data_pattern_map.find(dp);
	if(search == enum_to_str_data_pattern_map.end())
		return "UNKNOWN";
	else
		return search->second;
}

std::string einsim::get_all_possible_data_patterns(void)
{
	std::string ret;
	for(int s = einsim::DP_RANDOM; s < einsim::DP_UNKNOWN; s++)
		ret += einsim::enum_to_str_data_pattern((einsim::data_pattern)s) + ", ";
	ret.erase(ret.size() - 2, ret.size() - 1);
	return ret;
}

enum einsim::error_distribution einsim::str_to_enum_error_distribution(std::string str)
{
	static std::map< std::string, enum error_distribution > str_to_enum_error_distribution_map = 
	{
		  {"UNIFORM_RANDOM", ED_UNIFORM_RANDOM}
	};
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
	auto search = str_to_enum_error_distribution_map.find(str);
	if(search == str_to_enum_error_distribution_map.end())
		return ED_UNKNOWN;
	else
		return search->second;
}

std::string einsim::enum_to_str_error_distribution(enum einsim::error_distribution ed)
{
	static std::map< enum error_distribution, std::string > enum_to_str_error_distribution_map = 
	{
		  {ED_UNIFORM_RANDOM, "UNIFORM_RANDOM"}
	};
	auto search = enum_to_str_error_distribution_map.find(ed);
	if(search == enum_to_str_error_distribution_map.end())
		return "UNKNOWN";
	else
		return search->second;
}

std::string einsim::get_all_possible_error_distributions(void)
{
	std::string ret;
	for(int s = einsim::ED_UNIFORM_RANDOM; s < einsim::ED_UNKNOWN; s++)
		ret += einsim::enum_to_str_error_distribution((einsim::error_distribution)s) + ", ";
	ret.erase(ret.size() - 2, ret.size() - 1);
	return ret;
}

enum einsim::true_anti_cell_distribution einsim::str_to_enum_true_anti_cell_distribution(std::string str)
{
	static std::map< std::string, enum true_anti_cell_distribution > str_to_enum_true_anti_cell_distribution_map = 
	{
		  {"ALL_TRUE_OR_ALL_ANTI", CD_ALL_TRUE_OR_ALL_ANTI}
		, {"ALL_TRUE", CD_ALL_TRUE}
		, {"ALL_ANTI", CD_ALL_ANTI}
	};
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
	auto search = str_to_enum_true_anti_cell_distribution_map.find(str);
	if(search == str_to_enum_true_anti_cell_distribution_map.end())
		return CD_UNKNOWN;
	else
		return search->second;
}

std::string einsim::enum_to_str_true_anti_cell_distribution(enum einsim::true_anti_cell_distribution cd)
{
	static std::map< enum true_anti_cell_distribution, std::string > enum_to_str_true_anti_cell_distribution_map = 
	{
		  {CD_ALL_TRUE_OR_ALL_ANTI, "ALL_TRUE_OR_ALL_ANTI"}
		, {CD_ALL_TRUE, "ALL_TRUE"}
		, {CD_ALL_ANTI, "ALL_ANTI"}
	};
	auto search = enum_to_str_true_anti_cell_distribution_map.find(cd);
	if(search == enum_to_str_true_anti_cell_distribution_map.end())
		return "UNKNOWN";
	else
		return search->second;
}

std::string einsim::get_all_possible_true_anti_cell_distributions(void)
{
	std::string ret;
	for(int s = einsim::CD_ALL_TRUE_OR_ALL_ANTI; s < einsim::CD_UNKNOWN; s++)
		ret += einsim::enum_to_str_true_anti_cell_distribution((einsim::true_anti_cell_distribution)s) + ", ";
	ret.erase(ret.size() - 2, ret.size() - 1);
	return ret;
}

void einsim::inject(
	  Eigen::Matrix< ET, Eigen::Dynamic, 1 > &code_word
	, enum error_distribution ed
	, enum true_anti_cell_distribution cd
	, enum data_pattern dp
	, float rber
)
{
	if(!error_injector_initialized)
	{
		rng.seed(std::random_device{}());
		error_injector_initialized = true;
	}

	if(ed == ED_UNIFORM_RANDOM)
	{
		if(cd == CD_ALL_TRUE_OR_ALL_ANTI)
		{
			if(dp == DP_RANDOM)
			{
				std::bernoulli_distribution distribution(rber * 2.0); // 50% failure possibility -> all RBER must come from cells that can fail

				bool is_true_cell = ei_bernoulli_0_5(rng) ? 0 : 1;
				if(is_true_cell)
				{
					for(int i = 0; i < code_word.size(); i++)
						if(code_word(i))
							code_word(i) = distribution(rng) ? 0 : 1;
				}
				else
				{
					for(int i = 0; i < code_word.size(); i++)
						if(!code_word(i))
							code_word(i) = distribution(rng) ? 1 : 0;
				}
			}
			else if(dp == DP_ALL_ONES)
			{
				std::bernoulli_distribution distribution(rber * 2.0); // 50% failure possibility -> all RBER must come from cells that can fail

				bool is_true_cell = ei_bernoulli_0_5(rng) ? 0 : 1;
				if(is_true_cell)
					for(int i = 0; i < code_word.size(); i++)
						code_word(i) = distribution(rng) ? 0 : 1;
			}
			else if(dp == DP_CHARGED)
			{
				std::bernoulli_distribution distribution(rber * 2.0); // 50% failure possibility -> all RBER must come from cells that can fail

				// unfortunately, true-cells were determined at word_generation time.
				// fortunately, we can determine this by checking for all zeroes!
				bool is_true_cell = code_word.isZero(0) ? 0 : 1;
				if(is_true_cell)
				{
					for(int i = 0; i < code_word.size(); i++)
						if(code_word(i))
							code_word(i) = distribution(rng) ? 0 : 1;
				}
				else
				{
					for(int i = 0; i < code_word.size(); i++)
						if(!code_word(i))
							code_word(i) = distribution(rng) ? 1 : 0;
				}
			}
			else
			{
	            printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
	            exit(-1);					
			}
		}
		else
		{
            printf_both("[ERROR] Unsupported true-/anti-cell distribution: %d\n", (int)cd);
            exit(-1);
		}
	}
	else
	{
        printf_both("[ERROR] Unsupported error distribution: %d\n", (int)ed);
        exit(-1);
	}
}

// used only for testing/debugging, so we only handle the CHARGED and DISCHARGED cases
// and any cases that result in the same behavior (e.g., ALL_TRUE + ALL_ONES)
void einsim::inject_n(
	  Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
	, enum error_distribution ed
	, enum true_anti_cell_distribution cd
	, enum data_pattern dp
	, int n_errors
)
{
	if(!error_injector_initialized)
	{
		rng.seed(std::random_device{}());
		error_injector_initialized = true;
	}

	// Ideally, we need to create a mask of the bits that can actually fail
	// depending on the cd + dp. However, that's VERY slow. Instead, we only handle the cases
	// in which we effectively end up CHARGED or DISCHARGED
	bool is_charged = false;
	switch(cd)
	{
		case CD_ALL_TRUE_OR_ALL_ANTI:
		{
			switch(dp)
			{
				case DP_RANDOM:
					assert(0 && "there's no way to accurately compute this - the given parameters may or may not be impossible. use something deterministic or implement this under your own assumptions");
					break;

				case DP_ALL_ONES:
				{
					bool is_true_cell = ei_bernoulli_0_5(rng) ? 0 : 1;
					if(is_true_cell)
						is_charged = true;
					else
						is_charged = false;
					break;
				}

				case DP_CHARGED:
					is_charged = true;
					break;

				default:
					printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
					exit(-1);	
			}
			break;
		}

		case CD_ALL_TRUE:
		{
			switch(dp)
			{
				case DP_RANDOM:
					assert(0 && "there's no way to accurately compute this - the given parameters may or may not be impossible. use something deterministic or implement this under your own assumptions");
					break;

				case DP_ALL_ONES:
					is_charged = true;
					break;

				case DP_CHARGED:
					is_charged = true;
					break;

				default:
					printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
					exit(-1);	
			}
			break;
		}

		case CD_ALL_ANTI:
		{
			switch(dp)
			{
				case DP_RANDOM:
					assert(0 && "there's no way to accurately compute this - the given parameters may or may not be impossible. use something deterministic or implement this under your own assumptions");
					break;

				case DP_ALL_ONES:
					is_charged = false;
					break;

				case DP_CHARGED:
					is_charged = true;
					break;

				default:
					printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
					exit(-1);	
			}
			break;
		}

		default:
			printf_both("[ERROR] Unsupported true-/anti-cell distribution: %d\n", (int)cd);
			exit(-1);
			break;
	}

	// check if it's even possible to satisfy the request
	if(is_charged)
	{
		assert(n_errors <= word.size() && "impossible to meet requests for the given configuration");
		assert(ed == ED_UNIFORM_RANDOM && "unsupported error model requested!");

		// inject N errors uniform-randomly
		Eigen::Matrix< ET, Eigen::Dynamic, 1 > mask(word.size());
		int i;
		for(i = 0; i < n_errors; i++)
			mask[i] = 1;
		for(; i < word.size(); i++)
			mask[i] = 0;
		std::shuffle(mask.data(), mask.data() + mask.size(), rng_inject_errs);
		
		// inject the errors where requested
		for(int i = 0; i < word.size(); i++)
			word[i] ^= mask[i];
	}
	else
	{
		assert(n_errors == 0 && "impossible to meet requests for the given configuration");
		word = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(word.size());
		return;
	}
}