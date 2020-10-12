/**
 * @file observable.h
 *
 * @brief Module for measuring observables from simulation results
 *
 * This module is compartmentalized for logical reasons. However, the actual
 * measurements currently happen right as the simulation results are obtained in
 * the ``simulate'' and ``debug'' modules. If additional observables are introduced,
 * it will make sense to unify the measurement logic and move it here. 
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include <string>

namespace einsim
{
    /** @brief enumerates different observables that we can measure in simulation and experiment */ 
    enum observable
    {
          OBS_N_ERRORS_PER_BURST /**< counting the number of errors in a burst */
        , OBS_PER_BIT_ERROR_COUNT /**< histogram of how many times each bit fails */
        
        , OBS_UNKNOWN
    };

    /* routines for converting between the string and enum representations of error distributions */
    enum einsim::observable str_to_enum_observable(const std::string &str); /**< @brief converts observable string to enumeration */
    std::string enum_to_str_observable(enum einsim::observable obs); /**< @brief converts observable enumeration to a string */
    std::string get_all_possible_observables(void); /**< @brief returns a comma-separated string of all observable values */

}

#endif /* OBSERVABLE_H */