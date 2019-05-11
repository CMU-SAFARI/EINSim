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
        
        , OBS_MAX
    };

    /**
     * @brief enum to string conversion for the enum einsim::observable
     * 
     * @param obs enumeration value
     * @return string representation of the input value
     */
    std::string enum_to_str_observable(enum einsim::observable obs);
}

#endif /* OBSERVABLE_H */