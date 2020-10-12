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
#include <algorithm>

/* project includes */
#include "observable.h"

enum einsim::observable einsim::str_to_enum_observable(const std::string &str)
{
    static std::map< std::string, enum observable > str_to_enum_observable_map = 
    {
          {"N_ERRORS_PER_BURST", einsim::OBS_N_ERRORS_PER_BURST}
        , {"PER_BIT_ERROR_COUNT", einsim::OBS_PER_BIT_ERROR_COUNT}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_observable_map.find(str_uppercase);
    if(search == str_to_enum_observable_map.end())
        return OBS_UNKNOWN;
    else
        return search->second;
}

std::string einsim::enum_to_str_observable(enum einsim::observable obs)
{
    static std::map< enum einsim::observable, std::string > enum_to_str_observable_map = 
    {
          {einsim::OBS_N_ERRORS_PER_BURST, "N_ERRORS_PER_BURST"}
        , {einsim::OBS_PER_BIT_ERROR_COUNT, "PER_BIT_ERROR_COUNT"}
    };
    auto search = enum_to_str_observable_map.find(obs);
    if(search == enum_to_str_observable_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

std::string einsim::get_all_possible_observables(void)
{
    std::string ret;
    for(int s = einsim::OBS_N_ERRORS_PER_BURST; s < einsim::OBS_UNKNOWN; s++)
        ret += einsim::enum_to_str_observable((einsim::observable)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}
