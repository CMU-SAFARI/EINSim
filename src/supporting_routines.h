/**
 * @file supporting_routines.h
 *
 * @brief Miscellaneous utility variables and routines used throughout EINSim
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef SUPPORTING_ROUTINES_H
#define SUPPORTING_ROUTINES_H

#include <random>
#include <stdio.h>

/* libraries */
#include "Eigen/Eigen"

/** ET defines the Typename of all Eigen vectors/matrices used for ECC
  computations. 32-bit should be more than sufficient, but it's parameterized
  none the less. */
#define ET int32_t

/** convenience macro for performing a MOD2 on an entire Eigen vector */
#define MOD2_EIGEN_VEC(vector) ((vector).unaryExpr([](const ET x) { return x % 2; }))

extern int g_verbosity; /**< @brief globally controls the verbosity of the simulator */
extern FILE *g_output_file; /**< @brief FILE handle for the outfile so that all simulator parts can log to it */

/**
 * @brief computes the hamming distance between two binary Eigen vectors
 * 
 * @param a first input vector
 * @param b second input vector
 * @return hamming distance
 */
int hamming_distance(Eigen::Matrix< ET, Eigen::Dynamic, 1 > a, Eigen::Matrix< ET, Eigen::Dynamic, 1 > b);

/**
 * @brief prints to both stdout and the outfile given by the global output file
 * 
 * @param fmt format string
 * @return length of printed string
 */
int printf_both(const char *fmt, ...);

/**
 * @brief flushes both stdout and the outfile
 */
void fflush_both(void);

#endif /* SUPPORTING_ROUTINES_H	 */