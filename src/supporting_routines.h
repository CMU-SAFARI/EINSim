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
#include <initializer_list>
#include <set>
#include <iostream>
#include <stdio.h>

/* libraries */
#include "crc64/crc64.h"
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
 * @brief uses Gaussian elimination to row-reduce a matrix
 * 
 * @param m matrix to row-reduce
 * @param pivot_col pivot column to use, default 0
 * @return the matrix in RREF
 */
Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > row_reduce_to_rref(const Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &m, int pivot_col = 0);
void swap_matrix_rows(Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &m, int a, int b);


/**
 * @brief helper function to print ranges of values contained within a set
 * 
 * @param values set of numerical values to print
 * @tparam  T numerical type
 */
template< typename T >
void print_ranges(std::set< T > values)
{
    bool run_started = false;
    T last_val = 0;
    T min_e = *std::min_element(values.begin(), values.end());
    T max_e = *std::max_element(values.begin(), values.end());
    for(T i = min_e; i <= max_e + 1; i++)
    {
        if(values.find(i) == values.end())
        {
            if(run_started)
            {
                if(last_val == i - 1)
                    std::cout << last_val << " ";
                else
                    std::cout << last_val << "-" << i - 1 << " ";
                run_started = false;
            }
        }
        else if(!run_started)
        {
            last_val = i;
            run_started = true;
        }
    }
    std::cout << "]" << std::endl;
}

/**
 * @brief compute the crc64 hash function of a given list of matrices
 * 
 * @tparam T matrix type
 * @param mat_list list of matrices to hash
 * @return uint64_t hashed value
 */
template< typename T >
uint64_t hash_matrix(std::initializer_list< T > mat_list) 
{
    std::stringstream ss;
    for(const auto &mat : mat_list)
        for(int r = 0; r < mat.rows(); r++)
            for(int c = 0; c < mat.cols(); c++)
                ss << mat(r, c);
    return crc64(-1, ss.str().c_str(), ss.str().length());
};

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

/**
 * @brief convert 
 */

#endif /* SUPPORTING_ROUTINES_H  */