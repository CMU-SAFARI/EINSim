/**
 * @file bch_helpers.h
 *
 * @brief helper function prototypes for the BCH code implementation
 *     
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef BCH_HELPERS_H
#define BCH_HELPERS_H

#include "Eigen/Eigen"
#include "supporting_routines.h"

/**
 * @brief calculates the coefficients of the primitive polynomial specified by the function parameters
 * 
 * @param permutation selects between different pps that satisfy the gf_order requirements 
 * @param gf_order the order of the galois field for which to return the pp
 * @param primitive_polynomial coefficients of the primitive polynomial
 * @return 0 for success, !0 for error (e.g., no such polynomial available)
 */
int get_primitive_polynomial(int permutation, int gf_order, Eigen::Matrix< int, Eigen::Dynamic, 1 >& primitive_polynomial);

/**
 * @brief calculates the alpha_to and index_of tables for computing the GF per the primitive polynomial supplied
 * 
 * @param gf_order order of the galois field
 * @param primitive_polynomial coefficients of the primitive polynomial
 * @param alpha_to log table of GF(2**m) 
 * @param index_of anti-log table of GF(2**m) 
 */
void generate_gf(
      int gf_order
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > primitive_polynomial
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of);

/**
 * @brief attempts to calculate the generator polynomial for a BCH code of the desired parameters
 * 
 * @param m order of the galois field
 * @param code_len length of the BCH codeword
 * @param hd minimum hamming distance of the desired code
 * @param alpha_to log table of GF(2**m) 
 * @param index_of anti-log table of GF(2**m) 
 * @param g returns the coefficients of the BCH code generator polynomial
 * @param k returns the number of data bits supported by the BCH code
 * 
 * @return 0 on success, !0 if no such code is found
 */
int get_generator_polynomial(int m, int code_len, int hd 
    , const Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , const Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &g
    , int &k);

/**
 * @brief worker function to try to find an appropriate BCH code for the desired parameters
 * 
 * @param tid worker thread ID
 * @param permutation of all available primitive polynomials, which one to use
 * @param code_len length of the BCH codeword
 * @param nerrs_correctable how many errors the BCH code should be able to correct
 * @param gf_order order of the galois field
 * @param k returns the number of data bits supported by the BCH code
 * @param primitive_polynomial coefficients of the primitive polynomial
 * @param alpha_to log table of GF(2**m) 
 * @param index_of anti-log table of GF(2**m) 
 * @param generator_poly coefficients of the BCH code generator polynomial
 * 
 * @return 0 on success, !0 if no such code is found
 */
int get_bch_code_params(int tid, int permutation, int code_len, int nerrs_correctable, int gf_order
    , int &k
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &primitive_polynomial
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &generator_poly
);

/**
 * @brief wrapper function to try to find an appropriate BCH code for the desired parameters
 * 
 * @param permutation of all available primitive polynomials, which one to use
 * @param desired_n_data_bits desired length of the BCH codeword
 * @param nerrs_correctable desired correction capability of the BCH code
 * @param k returns the number of data bits supported by the BCH code
 * @param codeword_length returns the codeword length of the BCH code
 * @param gf_order returns the order of the galois field
 * @param primitive_polynomial coefficients of the primitive polynomial
 * @param alpha_to log table of GF(2**m) 
 * @param index_of anti-log table of GF(2**m) 
 * @param generator_poly coefficients of the BCH code generator polynomial
 * 
 * @return 0 on success, !0 if no such code is found
 */
int find_valid_bch_params(/*thread_pool &tp, */ int permutation, int desired_n_data_bits, int nerrs_correctable
    , int &k
    , int &codeword_length
    , int &gf_order
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &primitive_polynomial
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &generator_poly);

/**
 * @brief BCH encoding algorithm
 * 
 * @param length BCH codeword length 
 * @param k number of data bits supported by the BCH code
 * @param data_word_padded the dataword to be encoded, padded to the length of the codeword
 * @param g coefficients of the BCH code generator polynomial
 * 
 * @return encoded BCH codeword
 */
Eigen::Matrix< ET, Eigen::Dynamic, 1 > bch_encode(int length, int k
    , Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_padded
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > g);


/**
 * @brief BCH decoding algorithm
 * 
 * @param length BCH codeword length 
 * @param t correction capability of the code
 * @param n size of the multiplicative group of the GF order
 * @param code_word_padded the codeword to be decoded
 * @param alpha_to log table of GF(2**m) 
 * @param index_of anti-log table of GF(2**m) 
 * 
 * @return 0 on success, !0 if no such code is found
 */
void bch_decode(
      int length, int t, int n
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &code_word_padded
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &alpha_to
    , Eigen::Matrix< int, Eigen::Dynamic, 1 > &index_of);

#endif /* BCH_HELPERS_H */