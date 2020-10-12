/**
 * @file bch_helpers.cpp
 *
 * @brief Guts of the BCH encoding/decoding algorithms adapted from Robert
 *     Morelos-Zaragoza's bch3.h (original license below).
 *     
 * The original implementation has been dissected for use with EINSim,
 * including the following changes:
 *     - Primitive polynomials are taken from a hardcoded LUT instead of
 *       generated within a limited range. This enables having several, 
 *       easily configurable options and avoids unnecessary computation time
 *     - Internal data structures use Eigen data types instead of static
 *       storage arrays
 *     - Code in the original functions has been rearranged to fit the EINSim
 *       ecc_code class organization
 *     - C++ compatibility
 * 
 * Note that these functions should be threadsafe!
 *    
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
/*
 * File:    bch3.c
 * Title:   Encoder/decoder for binary BCH codes in C (Version 3.1)
 * Author:  Robert Morelos-Zaragoza
 * Date:    August 1994
 * Revised: June 13, 1997
 *
 * ===============  Encoder/Decoder for binary BCH codes in C =================
 *
 * Version 1:   Original program. The user provides the generator polynomial
 *              of the code (cumbersome!).
 * Version 2:   Computes the generator polynomial of the code.
 * Version 3:   No need to input the coefficients of a primitive polynomial of
 *              degree m, used to construct the Galois Field GF(2**m). The
 *              program now works for any binary BCH code of length such that:
 *              2**(m-1) - 1 < length <= 2**m - 1
 *
 * Note:        You may have to change the size of the arrays to make it work.
 *
 * The encoding and decoding methods used in this program are based on the
 * book "Error Control Coding: Fundamentals and Applications", by Lin and
 * Costello, Prentice Hall, 1983.
 *
 * Thanks to Patrick Boyle (pboyle@era.com) for his observation that 'bch2.c'
 * did not work for lengths other than 2**m-1 which led to this new version.
 * Portions of this program are from 'rs.c', a Reed-Solomon encoder/decoder
 * in C, written by Simon Rockliff (simon@augean.ua.oz.au) on 21/9/89. The
 * previous version of the BCH encoder/decoder in C, 'bch2.c', was written by
 * Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) on 5/19/92.
 *
 * NOTE:    
 *          The author is not responsible for any malfunctioning of
 *          this program, nor for any damage caused by it. Please include the
 *          original program along with these comments in any redistribution.
 *
 *  For more information, suggestions, or other ideas on implementing error
 *  correcting codes, please contact me at:
 *
 *                           Robert Morelos-Zaragoza
 *                           5120 Woodway, Suite 7036
 *                           Houston, Texas 77056
 *
 *                    email: r.morelos-zaragoza@ieee.org
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may 
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * == Copyright (c) 1994-7,  Robert Morelos-Zaragoza. All rights reserved.  ==
 *
 * m = order of the Galois field GF(2**m) 
 * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
 * length = length of the BCH code
 * t = error correcting capability (max. no. of errors the code corrects)
 * d = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
 * k = n - deg(g(x)) = dimension (no. of information bits/codeword) of the code
 * p[] = coefficients of a primitive polynomial used to generate GF(2**m)
 * g[] = coefficients of the generator polynomial, g(x)
 * alpha_to [] = log table of GF(2**m) 
 * index_of[] = antilog table of GF(2**m)
 * data[] = information bits = coefficients of data polynomial, i(x)
 * bb[] = coefficients of redundancy polynomial x^(length-k) i(x) modulo g(x)
 * numerr = number of errors 
 * errpos[] = error positions 
 * recd[] = coefficients of the received polynomial 
 * decerror = number of decoding errors (in _message_ positions) 
 *
 */
/* stdlib */
#include <iostream>
#include <array>
#include <string>
#include <set>
#include <vector>
#include <thread>

#include "supporting_routines.h"
#include "codes/bch_helpers.h"
#include "codes/bch_code.h"

static bool g_bch_helpers_print_errors = false;

/*
 *  gf_order (i.e., m): the degree of a primitive polynomial p(x) used to compute the
 *      Galois field GF(2**gf_order). Get precomputed coefficients p[] of p(x).
 */
int get_primitive_polynomial(int permutation, int gf_order, Eigen::Matrix< int, Eigen::Dynamic, 1 >& primitive_polynomial)
{
    if(gf_order < 3 || gf_order > 32)
        return -1;

    // primitive polynomials of the form (m, n, k ,t)
    // READ-only -> satisfies threadsafety
    static std::vector< std::vector< std::set< int > > > primitive_polynomial_entries = { 
          /* ( 3,    7,    4, 1) */ {{3, 1, 0}} // 0xb, o13
        , /* ( 4,   15,   11, 1) */ {{4, 1, 0}} // 0x13, o23
        , /* ( 5,   31,   26, 1) */ {{5, 2, 0}, {5, 4, 2, 1, 0}, {5, 4, 3, 2, 0}} // 0x25, o45
        , /* ( 6,   63,   57, 1) */ {{6, 1, 0}, {6, 5, 2, 1, 0}, {6, 5, 3, 2, 0}} // 0x83, o103
        , /* ( 7,  127,  120, 1) */ {{7, 1, 0}, {7, 3, 0}, {7, 3, 2, 1, 0}, {7, 4, 3, 2, 0}, {7, 5, 4, 3, 2, 1, 0}, {7, 6, 3, 1, 0}, {7, 6, 4, 2, 0}, {7, 6, 5, 2, 0}, {7, 6, 5, 4, 2, 1, 0}}
        , /* ( 8,  255,  247, 1) */ {{8, 4, 3, 2, 0}, {8, 5, 3, 1, 0}, {8, 6, 4, 3, 2, 1, 0}, {8, 6, 5, 1, 0}, {8, 6, 5, 2, 0}, {8, 6, 5, 3, 0}, {8, 7, 6, 1, 0}, {8, 7, 6, 5, 2, 1, 0}}
        , /* ( 9,  511,  502, 1) */ {{9, 4, 0}, {9, 5, 3, 2, 0}, {9, 6, 4, 3, 0}, {9, 6, 5, 3, 2, 1, 0}, {9, 6, 5, 4, 2, 1, 0}, {9, 7, 6, 4, 3, 1, 0}, {9, 8, 4, 1, 0}, {9, 8, 5, 4, 0}, {9, 8, 6, 5, 0}, {9, 8, 6, 5, 3, 1, 0}, {9, 8, 7, 2, 0}, {9, 8, 7, 3, 2, 1, 0}, {9, 8, 7, 6, 5, 1, 0}, {9, 8, 7, 6, 5, 3, 0}}
        , /* (10, 1023, 1013, 1) */ {{10, 3, 0}, {10, 4, 3, 1, 0}, {10, 6, 5, 3, 2, 1, 0}, {10, 8, 3, 2, 0}, {10, 8, 4, 3, 0}, {10, 8, 5, 1, 0}, {10, 8, 5, 4, 0}, {10, 8, 7, 6, 5, 2, 0}, {10, 8, 7, 6, 5, 4, 3, 1, 0}, {10, 9, 4, 1, 0}, {10, 9, 6, 5, 4, 3, 2, 1, 0}, {10, 9, 8, 6, 3, 2, 0}, {10, 9, 8, 6, 5, 1, 0}, {10, 9, 8, 7, 6, 5, 4, 3, 0}}
        , /* (  ,     ,     ,  ) */ {{11, 2, 0}, {11, 5, 3, 1, 0}, {11, 5, 3, 2, 0}, {11, 6, 5, 1, 0}, {11, 7, 3, 2, 0}, {11, 8, 5, 2, 0}, {11, 8, 6, 5, 4, 1, 0}, {11, 8, 6, 5, 4, 3, 2, 1, 0}, {11, 9, 4, 1, 0}, {11, 9, 8, 7, 4, 1, 0}, {11, 10, 3, 2, 0}, {11, 10, 7, 4, 3, 1, 0}, {11, 10, 8, 7, 5, 4, 3, 1, 0}, {11, 10, 9, 8, 3, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{12, 6, 4, 1, 0}, {12, 9, 3, 2, 0}, {12, 9, 8, 3, 2, 1, 0}, {12, 10, 9, 8, 6, 2, 0}, {12, 10, 9, 8, 6, 5, 4, 2, 0}, {12, 11, 6, 4, 2, 1, 0}, {12, 11, 9, 5, 3, 1, 0}, {12, 11, 9, 7, 6, 4, 0}, {12, 11, 9, 7, 6, 5, 0}, {12, 11, 9, 8, 7, 4, 0}, {12, 11, 9, 8, 7, 5, 2, 1, 0}, {12, 11, 10, 5, 2, 1, 0}, {12, 11, 10, 8, 6, 4, 3, 1, 0}, {12, 11, 10, 9, 8, 7, 5, 4, 3, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{13, 4, 3, 1, 0}, {13, 9, 7, 5, 4, 3, 2, 1, 0}, {13, 9, 8, 7, 5, 1, 0}, {13, 10, 9, 7, 5, 4, 0}, {13, 10, 9, 8, 6, 3, 2, 1, 0}, {13, 11, 8, 7, 4, 1, 0}, {13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}, {13, 12, 6, 5, 4, 3, 0}, {13, 12, 8, 7, 6, 5, 0}, {13, 12, 9, 8, 4, 2, 0}, {13, 12, 10, 8, 6, 4, 3, 2, 0}, {13, 12, 11, 5, 2, 1, 0}, {13, 12, 11, 8, 7, 6, 4, 1, 0}, {13, 12, 11, 9, 5, 3, 0}}
        , /* (  ,     ,     ,  ) */ {{14, 8, 6, 1, 0}, {14, 10, 6, 1, 0}, {14, 10, 9, 7, 6, 4, 3, 1, 0}, {14, 11, 6, 1, 0}, {14, 11, 9, 6, 5, 2, 0}, {14, 12, 9, 8, 7, 6, 5, 4, 0}, {14, 12, 11, 9, 8, 7, 6, 5, 3, 1, 0}, {14, 12, 11, 10, 9, 7, 4, 3, 0}, {14, 13, 6, 5, 3, 1, 0}, {14, 13, 10, 8, 7, 5, 4, 3, 2, 1, 0}, {14, 13, 11, 6, 5, 4, 2, 1, 0}, {14, 13, 11, 8, 5, 3, 2, 1, 0}, {14, 13, 12, 11, 10, 7, 6, 1, 0}, {14, 13, 12, 11, 10, 9, 6, 5, 0}}
        , /* (  ,     ,     ,  ) */ {{15, 1, 0}, {15, 4, 0}, {15, 7, 0}, {15, 7, 6, 3, 2, 1, 0}, {15, 10, 5, 1, 0}, {15, 10, 5, 4, 0}, {15, 10, 5, 4, 2, 1, 0}, {15, 10, 9, 7, 5, 3, 0}, {15, 10, 9, 8, 5, 3, 0}, {15, 11, 7, 6, 2, 1, 0}, {15, 12, 3, 1, 0}, {15, 12, 5, 4, 3, 2, 0}, {15, 12, 11, 8, 7, 6, 4, 2, 0}, {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 0}}
        , /* (  ,     ,     ,  ) */ {{16, 9, 8, 7, 6, 4, 3, 2, 0}, {16, 12, 3, 1, 0}, {16, 12, 7, 2, 0}, {16, 13, 12, 10, 9, 7, 6, 1, 0}, {16, 13, 12, 11, 7, 6, 3, 1, 0}, {16, 13, 12, 11, 10, 6, 2, 1, 0}, {16, 14, 10, 8, 3, 1, 0}, {16, 14, 13, 12, 6, 5, 3, 2, 0}, {16, 14, 13, 12, 10, 7, 0}, {16, 15, 10, 6, 5, 3, 2, 1, 0}, {16, 15, 11, 9, 8, 7, 5, 4, 2, 1, 0}, {16, 15, 11, 10, 7, 6, 5, 3, 2, 1, 0}, {16, 15, 11, 10, 9, 6, 2, 1, 0}, {16, 15, 11, 10, 9, 8, 6, 4, 2, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{17, 3, 0}, {17, 3, 2, 1, 0}, {17, 5, 0}, {17, 6, 0}, {17, 8, 4, 3, 0}, {17, 8, 7, 6, 4, 3, 0}, {17, 10, 9, 8, 6, 5, 3, 2, 0}, {17, 12, 6, 3, 2, 1, 0}, {17, 12, 9, 5, 4, 3, 2, 1, 0}, {17, 12, 9, 7, 6, 4, 3, 2, 0}, {17, 14, 11, 7, 5, 3, 2, 1, 0}, {17, 15, 13, 11, 9, 7, 5, 3, 0}, {17, 15, 13, 11, 9, 7, 6, 4, 2, 1, 0}, {17, 16, 3, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{18, 5, 4, 3, 2, 1, 0}, {18, 7, 0}, {18, 7, 5, 2, 1, 0}, {18, 8, 2, 1, 0}, {18, 9, 7, 6, 5, 4, 0}, {18, 9, 8, 6, 5, 4, 2, 1, 0}, {18, 9, 8, 7, 6, 4, 2, 1, 0}, {18, 10, 7, 5, 0}, {18, 10, 8, 5, 0}, {18, 10, 8, 7, 6, 5, 4, 3, 2, 1, 0}, {18, 10, 9, 3, 0}, {18, 13, 6, 4, 0}, {18, 15, 5, 2, 0}, {18, 15, 9, 2, 0}}
        , /* (  ,     ,     ,  ) */ {{19, 5, 2, 1, 0}, {19, 5, 4, 3, 2, 1, 0}, {19, 6, 2, 1, 0}, {19, 6, 5, 3, 2, 1, 0}, {19, 6, 5, 4, 3, 2, 0}, {19, 7, 5, 3, 2, 1, 0}, {19, 8, 7, 5, 0}, {19, 8, 7, 5, 4, 3, 2, 1, 0}, {19, 8, 7, 6, 4, 3, 2, 1, 0}, {19, 9, 8, 5, 0}, {19, 9, 8, 6, 5, 3, 2, 1, 0}, {19, 9, 8, 7, 4, 3, 2, 1, 0}, {19, 11, 9, 8, 7, 6, 5, 4, 3, 2, 0}, {19, 11, 10, 8, 7, 5, 4, 3, 2, 1, 0}, {19, 16, 13, 10, 7, 4, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{20, 3, 0}, {20, 9, 5, 3, 0}, {20, 11, 8, 6, 3, 2, 0}, {20, 14, 10, 9, 8, 6, 5, 4, 0}, {20, 17, 14, 10, 7, 4, 3, 2, 0}, {20, 19, 4, 3, 0}}
        , /* (  ,     ,     ,  ) */ {{21, 2, 0}, {21, 8, 7, 4, 3, 2, 0}, {21, 10, 6, 4, 3, 2, 0}, {21, 13, 5, 2, 0}, {21, 14, 7, 2, 0}, {21, 14, 7, 6, 3, 2, 0}, {21, 14, 12, 7, 6, 4, 3, 2, 0}, {21, 15, 10, 9, 5, 4, 3, 2, 0}, {21, 20, 19, 18, 5, 4, 3, 2, 0}}
        , /* (  ,     ,     ,  ) */ {{22, 1, 0}, {22, 9, 5, 1, 0}, {22, 14, 13, 12, 7, 3, 2, 1, 0}, {22, 17, 9, 7, 2, 1, 0}, {22, 17, 13, 12, 8, 7, 2, 1, 0}, {22, 20, 18, 16, 6, 4, 2, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{23, 5, 0}, {23, 5, 4, 1, 0}, {23, 11, 10, 7, 6, 5, 0}, {23, 12, 5, 4, 0}, {23, 15, 10, 9, 7, 5, 4, 3, 0}, {23, 16, 13, 6, 5, 3, 0}, {23, 17, 11, 5, 0}, {23, 17, 11, 9, 8, 5, 4, 1, 0}, {23, 18, 16, 13, 11, 8, 5, 2, 0}, {23, 21, 7, 5, 0}}
        , /* (  ,     ,     ,  ) */ {{24, 7, 2, 1, 0}, {24, 21, 19, 18, 17, 16, 15, 14, 13, 10, 9, 5, 4, 1, 0}, {24, 22, 20, 18, 16, 14, 11, 9, 8, 7, 5, 4, 0}}
        , /* (  ,     ,     ,  ) */ {{25, 3, 0}, {25, 3, 2, 1, 0}, {25, 11, 9, 8, 6, 4, 3, 2, 0}, {25, 12, 4, 3, 0}, {25, 12, 11, 8, 7, 6, 4, 3, 0}, {25, 17, 10, 3, 2, 1, 0}, {25, 18, 12, 11, 6, 5, 4, 3, 0}, {25, 20, 5, 3, 0}, {25, 20, 16, 11, 5, 3, 2, 1, 0}, {25, 23, 21, 19, 9, 7, 5, 3, 0}}
        , /* (  ,     ,     ,  ) */ {{26, 6, 2, 1, 0}, {26, 19, 16, 15, 14, 13, 11, 9, 8, 7, 6, 5, 3, 2, 0}, {26, 21, 18, 16, 15, 13, 12, 11, 9, 8, 6, 5, 4, 3, 0}, {26, 22, 20, 19, 16, 13, 11, 9, 8, 7, 5, 4, 2, 1, 0}, {26, 22, 21, 16, 12, 11, 10, 8, 5, 4, 3, 1, 0}, {26, 23, 22, 21, 19, 18, 15, 14, 13, 11, 10, 9, 8, 6, 5, 2, 0}, {26, 24, 21, 17, 16, 14, 13, 11, 7, 6, 4, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{27, 5, 2, 1, 0}, {27, 18, 11, 10, 9, 5, 4, 3, 0}, {27, 22, 13, 11, 6, 5, 4, 3, 0}, {27, 22, 17, 15, 14, 13, 6, 1, 0}, {27, 22, 21, 20, 18, 17, 15, 13, 12, 7, 5, 0}, {27, 24, 19, 16, 12, 8, 7, 3, 2, 1, 0}, {27, 24, 21, 19, 16, 13, 11, 9, 6, 5, 4, 3, 0}, {27, 25, 23, 21, 13, 11, 9, 8, 7, 6, 5, 3, 2, 1, 0}, {27, 25, 23, 21, 20, 19, 18, 16, 14, 10, 8, 7, 4, 3, 0}}
        , /* (  ,     ,     ,  ) */ {{28, 3, 0}, {28, 13, 11, 9, 5, 3, 0}, {28, 18, 17, 16, 9, 5, 4, 3, 0}, {28, 19, 17, 15, 10, 6, 3, 2, 0}, {28, 22, 11, 10, 4, 3, 0}, {28, 24, 20, 16, 12, 8, 4, 3, 0}}
        , /* (  ,     ,     ,  ) */ {{29, 2, 0}, {29, 12, 7, 2, 0}, {29, 18, 14, 6, 3, 2, 0}, {29, 19, 16, 6, 3, 2, 0}, {29, 20, 11, 2, 0}, {29, 20, 16, 11, 8, 4, 3, 2, 0}, {29, 21, 5, 2, 0}, {29, 23, 10, 9, 5, 4, 3, 2, 0}, {29, 24, 14, 13, 8, 4, 3, 2, 0}, {29, 26, 5, 2, 0}}
        , /* (  ,     ,     ,  ) */ {{30, 23, 2, 1, 0}, {30, 24, 20, 16, 14, 13, 11, 7, 2, 1, 0}, {30, 24, 21, 20, 18, 15, 13, 12, 9, 7, 6, 4, 3, 1, 0}, {30, 25, 24, 23, 19, 18, 16, 14, 11, 8, 6, 4, 3, 1, 0}, {30, 27, 25, 24, 23, 22, 19, 16, 12, 10, 8, 7, 6, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{31, 3, 0}, {31, 3, 2, 1, 0}, {31, 13, 8, 3, 0}, {31, 16, 8, 4, 3, 2, 0}, {31, 20, 15, 5, 4, 3, 0}, {31, 20, 18, 7, 5, 3, 0}, {31, 21, 12, 3, 2, 1, 0}, {31, 23, 22, 15, 14, 7, 4, 3, 0}, {31, 25, 19, 14, 7, 3, 2, 1, 0}, {31, 27, 23, 19, 15, 11, 7, 3, 0}, {31, 27, 23, 19, 15, 11, 10, 9, 7, 6, 5, 3, 2, 1, 0}}
        , /* (  ,     ,     ,  ) */ {{32, 22, 2, 1, 0}, {32, 22, 21, 20, 18, 17, 15, 13, 12, 10, 8, 6, 4, 1, 0}, {32, 23, 17, 16, 14, 10, 8, 7, 6, 5, 3, 0}, {32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1, 0}, {32, 27, 26, 25, 24, 23, 22, 17, 13, 11, 10, 9, 8, 7, 2, 1, 0}, {32, 28, 19, 18, 16, 14, 11, 10, 9, 6, 5, 1, 0}}
    };

    primitive_polynomial = Eigen::Matrix< int, Eigen::Dynamic, 1 >::Zero(gf_order + 1);
    std::vector< std::set< int > > &elements_array = primitive_polynomial_entries[gf_order - 3];
    std::set< int > &elements = elements_array[permutation % elements_array.size()];
    for(const auto& element : elements)
        primitive_polynomial[element] = 1;
    return 0;
}


/*
 * Generate field GF(2**m) from the irreducible polynomial p(X) with
 * coefficients in p[0]..p[m].
 *
 * Lookup tables:
 *   index->polynomial form: alpha_to[] contains j=alpha^i;
 *   polynomial form -> index form: index_of[j=alpha^i] = i
 *
 * alpha=2 is the primitive element of GF(2**m) 
 */
void generate_gf(
      int gf_order
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > primitive_polynomial
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of)
{
    int n = (1 << gf_order) - 1; // size of multiplicative group GF(2^m)
    int    i, mask;
    alpha_to = Eigen::Matrix< int, Eigen::Dynamic, 1 >::Zero(n);
    index_of = Eigen::Matrix< int, Eigen::Dynamic, 1 >::Zero(n + 1);

    mask = 1;
    alpha_to[gf_order] = 0;
    for (i = 0; i < gf_order; i++) {
        alpha_to[i] = mask;
        index_of[alpha_to[i]] = i;
        if (primitive_polynomial[i] != 0)
            alpha_to[gf_order] ^= mask;
        mask <<= 1;
    }
    index_of[alpha_to[gf_order]] = gf_order;
    mask >>= 1;
    for (i = gf_order + 1; i < n; i++) {
        if (alpha_to[i - 1] >= mask)
          alpha_to[i] = alpha_to[gf_order] ^ ((alpha_to[i - 1] ^ mask) << 1);
        else
          alpha_to[i] = alpha_to[i - 1] << 1;
        index_of[alpha_to[i]] = i;
        // printf("n: %d MAX: %d\n", n, alpha_to[i]);
    }
    index_of[0] = -1;
}

/*
 * Compute the generator polynomial of a binary BCH code. Fist generate the
 * cycle sets modulo 2**m - 1, cycle[][] = (i, 2*i, 4*i, ..., 2^l*i). Then
 * determine those cycle sets that contain integers in the set of (d-1)
 * consecutive integers {1..(d-1)}. The generator polynomial is calculated
 * as the product of linear factors of the form (x+alpha^i), for every i in
 * the above cycle sets.
 */
int get_generator_polynomial(int m, int code_len, int hd 
    , const Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , const Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &g
    , int &k)
{
    int n = (1 << m) - 1; // size of the multiplicative group of GF(2**m)
    int    ii, jj, ll, kaux;
    int    test, aux, nocycles, root, noterms, rdncy;
    Eigen::Matrix< int, Eigen::Dynamic, 1 > size;
    Eigen::Matrix< int, Eigen::Dynamic, 1 > min;
    Eigen::Matrix< int, Eigen::Dynamic, 1 > zeros;
    Eigen::Matrix< int, Eigen::Dynamic, 32 > cycle;
    size.resize(n, 1);
    min.resize(n, 1);
    zeros.resize(n, 1);
    cycle.resize(n, 32);

    /* Generate cycle sets modulo n (n = 2**m - 1) */
    cycle(0, 0) = 0;
    size[0] = 1;
    cycle(1, 0) = 1;
    size[1] = 1;
    jj = 1;         /* cycle set index */
    if (m > 9) 
        printf("[WARNING] Computing cycle sets modulo %d - this is slow\n", n);

    do {
        /* Generate the jj-th cycle set */
        ii = 0;
        do {
            ii++;
            cycle(jj, ii) = (cycle(jj, ii - 1) * 2) % n;
            size[jj]++;
            aux = (cycle(jj, ii) * 2) % n;
        } while (aux != cycle(jj, 0));
        /* Next cycle set representative */
        ll = 0;
        do {
            ll++;
            test = 0;
            for (ii = 1; ((ii <= jj) && (!test)); ii++) 
            /* Examine previous cycle sets */
              for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
                 if (ll == cycle(ii, kaux))
                    test = 1;
        } while ((test) && (ll < (n - 1)));
        if (!(test)) {
            jj++;   /* next cycle set index */
            cycle(jj, 0) = ll;
            size[jj] = 1;
        }
    } while (ll < (n - 1));
    nocycles = jj;      /* number of cycle sets modulo n */

    /* Search for roots 1, 2, ..., d-1 in cycle sets */
    kaux = 0;
    rdncy = 0;
    for (ii = 1; ii <= nocycles; ii++)
    {
        min[kaux] = 0;
        test = 0;
        for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
            for (root = 1; ((root < hd) && (!test)); root++)
                if (root == cycle(ii, jj))
                {
                    test = 1;
                    min[kaux] = ii;
                }
        if (min[kaux])
        {
            rdncy += size[min[kaux]];
            kaux++;
        }
    }
    noterms = kaux;
    kaux = 1;
    for (ii = 0; ii < noterms; ii++)
        for (jj = 0; jj < size[min[ii]]; jj++)
        {
            zeros[kaux] = cycle(min[ii], jj);
            kaux++;
        }

    k = code_len - rdncy;

    // if we have no (or fewer) data bits, this is an invalid code
    if(k <= 0)
    {
        if(g_bch_helpers_print_errors)
            printf("ERROR: Invalid parameters: n: %d, k: %d!\n", code_len, k);
        return -1;
    }

    // printf("This is a (n: %d, k: %d, d: %d) binary BCH code\n", code_len, k, hd);

    /* Compute the generator polynomial */
    g = Eigen::Matrix< int, Eigen::Dynamic, 1 >::Zero(std::max(rdncy + 1, 2));
    g[0] = alpha_to[zeros[1]];
    g[1] = 1;       /* g(x) = (X + zeros[1]) initially */
    for(ii = 2; ii <= rdncy; ii++)
    {
        g[ii] = 1;
        for(jj = ii - 1; jj > 0; jj--)
            if(g[jj] != 0)
                g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % n];
            else
                g[jj] = g[jj - 1];
        g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % n];
    }

    // printf("Generator polynomial:\ng(x) = ");
    // for (ii = 0; ii <= rdncy; ii++)
    // {
    //     printf("%d", g[ii]);
    //     if (ii && ((ii % 50) == 0))
    //         printf("\n");
    // }
    // printf("\n");

    return 0;
}

// returns the number of data bits in such a code
int get_bch_code_params(int tid, int permutation, int code_len, int nerrs_correctable, int gf_order
    , int &k
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &primitive_polynomial
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &generator_poly
)
{
    int hamming_dist = 2 * nerrs_correctable + 1;

    /* sanity-check that t < 2^(m - 1) */
    if(nerrs_correctable >= (1 << (gf_order - 1)))
    {
        printf("Rejecting invalid BCH code of perm: %d m: %d t: %d\n", permutation, gf_order, nerrs_correctable);
        return -1;                
    }

    /* get the primitive polynomial */
    if(get_primitive_polynomial(permutation, gf_order, primitive_polynomial))
    {
        printf("Rejecting invalid BCH code of perm: %d m: %d n: %d, t: %d\n", permutation, gf_order, code_len, nerrs_correctable);
        return -1;
    }

    /* generate the entire GF associated with this primitive polynomial */
    generate_gf(gf_order, primitive_polynomial, alpha_to, index_of);

    /* construct the generator polynomial */
    if(get_generator_polynomial(gf_order, code_len, hamming_dist, alpha_to, index_of, generator_poly, k))
    {
        printf("Rejecting invalid BCH code of perm: %d m: %d n: %d, k: %d t: %d\n", permutation, gf_order, code_len, k, nerrs_correctable);
        return -1;
    }

    if(g_verbosity > 0)
        printf("Considering valid BCH code of perm: %d m: %d n: %d, k: %d t: %d\n", permutation, gf_order, code_len, k, nerrs_correctable);
    return 0;
}

int find_valid_bch_params(/*thread_pool &tp, */ int permutation, int desired_n_data_bits, int nerrs_correctable
    , int &k
    , int &codeword_length
    , int &gf_order
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &primitive_polynomial
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &generator_poly)
{
    // find the requried GF order to support this nubmer of bits
    int min_m = 1;
    while((1 << min_m) < desired_n_data_bits)  // implements Po2 ceil
        min_m++;

    // find the minimum gf_order that will work for this case
    for(int m = min_m; m <= 13; m++)
    {
        int n = (1 << m) - 1;

        if(get_bch_code_params(1, permutation, n, nerrs_correctable, m
            , k, primitive_polynomial, alpha_to, index_of, generator_poly))
            continue;
        if(k >= desired_n_data_bits)
        {
            codeword_length = n;
            gf_order = m;
            return 0;
        }
    }
    return -1;
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > bch_encode(int length, int k
    , Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_padded
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > g)
{
    /*
     * Compute redundacy bb[], the coefficients of b(x). The redundancy
     * polynomial b(x) is the remainder after dividing x^(length-k)*data(x)
     * by the generator polynomial g(x).
     */
    int    i, j;
    int    feedback;
    
    /*
     * rpoly[] are the coefficients of the redundancy polynomial c(x) = x**(length-k)*data(x) + b(x)
     */
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > rpoly = Eigen::Matrix< ET, Eigen::Dynamic, 1 >::Zero(length - k);
    for(i = k - 1; i >= 0; i--)
    {
        feedback = data_word_padded[i] ^ rpoly[length - k - 1];
        if(feedback != 0)
        {
            for(j = length - k - 1; j > 0; j--)
                if(g[j] != 0)
                    rpoly[j] = rpoly[j - 1] ^ feedback;
                else
                    rpoly[j] = rpoly[j - 1];
            rpoly[0] = g[0] && feedback;
        }
        else
        {
            for(j = length - k - 1; j > 0; j--)
                rpoly[j] = rpoly[j - 1];
            rpoly[0] = 0;
        }
    }

    // for(i = 0; i < k; i++)
    //     rpoly[i + length - k] = data_word_padded[i];
    return rpoly;
}

/** 
 * @brief local variables necessary for bch_decode
 * 
 * these are too big to allocate on an OSX default pthread 512KiB stack :(
 * therefore, we pull them out into a single POD struct and dynamically
 * allocate it to ensure that calls to bch_decode are threadsafe
 */
struct bch_decode_locals
{
    int elp[1024][1024];
    int d[1026];
    int l[1026];
    int u_lu[1026]; 
    int s[1025];
    int loc[200];
    int reg[201];
};

/*
 * Simon Rockliff's implementation of Berlekamp's algorithm.
 *
 * Assume we have received bits in code_word[i], i=0..(n-1).
 *
 * Compute the 2*t syndromes by substituting alpha^i into rec(X) and
 * evaluating, storing the syndromes in s[i], i=1..2t (leave s[0] zero) .
 * Then we use the Berlekamp algorithm to find the error location polynomial
 * elp[i].
 *
 * If the degree of the elp is >t, then we cannot correct all the errors, and
 * we have detected an uncorrectable error pattern. We output the information
 * bits uncorrected.
 *
 * If the degree of elp is <=t, we substitute alpha^i , i=1..n into the elp
 * to get the roots, hence the inverse roots, the error location numbers.
 * This step is usually called "Chien's search".
 *
 * If the number of errors located is not equal the degree of the elp, then
 * the decoder assumes that there are more than t errors and cannot correct
 * them, only detect them. We output the information bits uncorrected.
 */
void bch_decode(
      int length, int t, int n
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &code_word_padded
    , const Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , const Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of)
{
    int    i, j, u, q, t2, count = 0, syn_error = 0;
    struct bch_decode_locals *bdl = new struct bch_decode_locals;
    
    // std::thread::id this_id = std::this_thread::get_id();
    // std::cout << "my tid: " << this_id << std::endl;
    // std::map< std::thread::id, 

    t2 = 2 * t;

    /* first form the syndromes */
    // printf("S(x) = ");
    for(i = 1; i <= t2; i++)
    {
        bdl->s[i] = 0;
        for(j = 0; j < length; j++)
            if(code_word_padded[j] != 0)
                bdl->s[i] ^= alpha_to[(i * j) % n];
        if(bdl->s[i] != 0)
            syn_error = 1; /* set error flag if non-zero syndrome */
        /*
         * Note:    If the code is used only for ERROR DETECTION, then
         *          exit program here indicating the presence of errors.
         */
        /* convert syndrome from polynomial form to index form  */
        bdl->s[i] = index_of[bdl->s[i]];
        // printf("%3d ", s[i]);
    }
    // printf("\n");

    if(syn_error)       /* if there are errors, try to correct them */
    {
        /*
         * Compute the error location polynomial via the Berlekamp
         * iterative algorithm. Following the terminology of Lin and
         * Costello's book :   d[u] is the 'mu'th discrepancy, where
         * u='mu'+1 and 'mu' (the Greek letter!) is the step number
         * ranging from -1 to 2*t (see L&C),  l[u] is the degree of
         * the elp at that step, and u_l[u] is the difference between
         * the step number and the degree of the elp.
         */
        /* initialise table entries */
        bdl->d[0] = 0;           /* index form */
        bdl->d[1] = bdl->s[1];        /* index form */
        bdl->elp[0][0] = 0;      /* index form */
        bdl->elp[1][0] = 1;      /* polynomial form */
        for(i = 1; i < t2; i++)
        {
            bdl->elp[0][i] = -1; /* index form */
            bdl->elp[1][i] = 0;  /* polynomial form */
        }
        bdl->l[0] = 0;
        bdl->l[1] = 0;
        bdl->u_lu[0] = -1;
        bdl->u_lu[1] = 0;
        u = 0;

        do
        {
            u++;
            if(bdl->d[u] == -1)
            {
                bdl->l[u + 1] = bdl->l[u];
                for(i = 0; i <= bdl->l[u]; i++)
                {
                    bdl->elp[u + 1][i] = bdl->elp[u][i];
                    bdl->elp[u][i] = index_of[bdl->elp[u][i]];
                }
            }
            else
                /*
                 * search for words with greatest u_lu[q] for
                 * which d[q]!=0
                 */
            {
                q = u - 1;
                while((bdl->d[q] == -1) && (q > 0))
                    q--;
                /* have found first non-zero d[q]  */
                if(q > 0)
                {
                    j = q;
                    do
                    {
                        j--;
                        if((bdl->d[j] != -1) && (bdl->u_lu[q] < bdl->u_lu[j]))
                            q = j;
                    }
                    while(j > 0);
                }

                /*
                 * have now found q such that d[u]!=0 and
                 * u_lu[q] is maximum
                 */
                /* store degree of new elp polynomial */
                if(bdl->l[u] > bdl->l[q] + u - q)
                    bdl->l[u + 1] = bdl->l[u];
                else
                    bdl->l[u + 1] = bdl->l[q] + u - q;

                /* form new elp(x) */
                for(i = 0; i < t2; i++)
                    bdl->elp[u + 1][i] = 0;
                for(i = 0; i <= bdl->l[q]; i++)
                    if(bdl->elp[q][i] != -1)
                        bdl->elp[u + 1][i + u - q] =
                            alpha_to[(bdl->d[u] + n - bdl->d[q] + bdl->elp[q][i]) % n];
                for(i = 0; i <= bdl->l[u]; i++)
                {
                    bdl->elp[u + 1][i] ^= bdl->elp[u][i];
                    bdl->elp[u][i] = index_of[bdl->elp[u][i]];
                }
            }
            bdl->u_lu[u + 1] = u - bdl->l[u + 1];

            /* form (u+1)th discrepancy */
            if(u < t2)
            {
                /* no discrepancy computed on last iteration */
                if(bdl->s[u + 1] != -1)
                    bdl->d[u + 1] = alpha_to[bdl->s[u + 1]];
                else
                    bdl->d[u + 1] = 0;
                for(i = 1; i <= bdl->l[u + 1]; i++)
                    if((bdl->s[u + 1 - i] != -1) && (bdl->elp[u + 1][i] != 0))
                        bdl->d[u + 1] ^= alpha_to[(bdl->s[u + 1 - i]
                                              + index_of[bdl->elp[u + 1][i]]) % n];
                /* put d[u+1] into index form */
                bdl->d[u + 1] = index_of[bdl->d[u + 1]];
            }
        }
        while((u < t2) && (bdl->l[u + 1] <= t));

        u++;
        if(bdl->l[u] <= t)   /* Can correct errors */
        {
            /* put elp into index form */
            for(i = 0; i <= bdl->l[u]; i++)
                bdl->elp[u][i] = index_of[bdl->elp[u][i]];

            // printf("sigma(x) = ");
            // for(i = 0; i <= l[u]; i++)
            //     printf("%3d ", bdl->elp[u][i]);
            // printf("\n");
            // printf("Roots: ");

            /* Chien search: find roots of the error location polynomial */
            for(i = 1; i <= bdl->l[u]; i++)
                bdl->reg[i] = bdl->elp[u][i];
            count = 0;
            for(i = 1; i <= n; i++)
            {
                q = 1;
                for(j = 1; j <= bdl->l[u]; j++)
                    if(bdl->reg[j] != -1)
                    {
                        bdl->reg[j] = (bdl->reg[j] + j) % n;
                        q ^= alpha_to[bdl->reg[j]];
                    }
                if(!q)
                {
                    /* store root and error
                             * location number indices */
                    // root[count] = i;
                    bdl->loc[count] = n - i;
                    count++;
            //        printf("%3d ", n - i);
                }
            }
            // printf("\n");
            if(count == bdl->l[u])
                /* no. roots = degree of elp hence <= t errors */
                for(i = 0; i < bdl->l[u]; i++)
                {
            //        std::cout << "Loc[" << i << "] = " << bdl->loc[i] << std::endl;
                    code_word_padded[bdl->loc[i]] ^= 1;
                }
            // else    /* elp has degree >t hence cannot solve */
            //     printf("Incomplete decoding: errors detected\n");
        }
    }

    delete bdl;
}
