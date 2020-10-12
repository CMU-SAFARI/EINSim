/**
 * @file supporting_routines.cpp
 *
 * @brief Miscellaneous helper routines
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <random>
#include <iostream>
#include <stdarg.h>

/* libraries */
#include "Eigen/Eigen"

/* project includes */
#include "supporting_routines.h"

int g_verbosity = 0;
FILE *g_output_file = NULL;

int printf_both(const char *fmt, ...)
{
    int ret = 0;
    {
        va_list ap;
        va_start(ap, fmt);
        ret = vprintf(fmt, ap);
        va_end(ap);
    }
    if(g_output_file != stdout)
    {
        va_list ap;
        va_start(ap, fmt);
        ret = vfprintf(g_output_file, fmt, ap);   
        va_end(ap);
    }
    return ret;
}

void fflush_both(void)
{
    fflush(stdout);
    assert(g_output_file);
    fflush(g_output_file);
}

/**
 * @brief compute the Hamming distance between two vectors
 * 
 * @param a vector a
 * @param b vector b
 * @return int Hammign distance
 */
int hamming_distance(Eigen::Matrix< ET, Eigen::Dynamic, 1 > a, Eigen::Matrix< ET, Eigen::Dynamic, 1 > b)
{
    return (a - b).count();
}

/**
 * @brief swap two matrix rows
 * 
 * @param m matrix
 * @param a index of row a
 * @param b index of row b
 */
void swap_matrix_rows(Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &m, int a, int b)
{
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > temp = m.row(a);
    m.row(a) = m.row(b);
    m.row(b) = temp;
}

/**
 * @brief add a row to another row of a matrix
 * 
 * @param m matrix
 * @param a source row index
 * @param b destination row index
 */
void add_row_a_to_row_b(Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &m, int a, int b)
{
    m.row(b) = MOD2_EIGEN_VEC(m.row(a) + m.row(b));
}


Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > row_reduce_to_rref(const Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &m, int pivot_col)
{
    Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > ret(m);

    int nrows = m.rows();
    int ncols = m.cols();

    int cur_pivot = 0;
    for(int c = pivot_col; c < ncols; c++)
    {
        // std::cout << "[RREF] column: " << c << std::endl << ret << std::endl;

        // find first nonzero row in this column
        for(int r = cur_pivot; r < nrows; r++)
        {
            if(ret(r, c) == 1)
            {
                // swap the nonzero row to the pivot position if necessary
                if(r != cur_pivot)
                {
                    // std::cout << "[RREF] swapping rows " << r << " and " << cur_pivot << std::endl;
                    swap_matrix_rows(ret, r, cur_pivot);
                }

                // rearrange/add rows such that the r is the only nonzero in this column
                for(int r = 0; r < nrows; r++)
                {
                    if(r != cur_pivot && ret(r, c) == 1)
                    {
                        // std::cout << "[RREF] adding row " << cur_pivot << " to row " << r << std::endl;
                        add_row_a_to_row_b(ret, cur_pivot, r);
                    }
                }
                cur_pivot++;
                break;
            }
        }
    }

    return ret;
}