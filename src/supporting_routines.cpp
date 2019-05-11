/**
 * @file supporting_routines.cpp
 *
 * @brief Miscellaneous helper routines
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <random>
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

int hamming_distance(Eigen::Matrix< ET, Eigen::Dynamic, 1 > a, Eigen::Matrix< ET, Eigen::Dynamic, 1 > b)
{
	return (a - b).count();
}