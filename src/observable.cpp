/**
 * @file observable.cpp
 *
 * @brief Definitions for measuring observables from simulation results
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <string>
#include <set>
#include <map>

/* project includes */
#include "observable.h"

std::string einsim::enum_to_str_observable(enum einsim::observable obs)
{
	static std::map< enum einsim::observable, std::string > enum_to_str_observable_map = 
	{
		  {einsim::OBS_N_ERRORS_PER_BURST, "N_ERRORS_PER_BURST"}
	};
	auto search = enum_to_str_observable_map.find(obs);
	if(search == enum_to_str_observable_map.end())
		return "UNKNOWN";
	else
		return search->second;
}