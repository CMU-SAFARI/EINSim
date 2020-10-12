/**
 * @file bch_code.h
 *
 * @brief Binary BCH code implementation
 * 
 * This implementation is based largely on the BCH encoding/decoding
 * algorithms adapted from Robert Morelos-Zaragoza's bch3.h.
 * ``bch_helpers.cpp'' contains the implementation and original license for
 * the adapted code.
 *     
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef BCH_CODE_H
#define BCH_CODE_H

#include "Eigen/Eigen"
#include "supporting_routines.h"
#include "codes/bch_helpers.h"
#include "ecc_code.h"

namespace einsim
{
    /**
     * @brief implements a binary BCH code
     * 
     * The BCH code is specified using:
     *      * int n_errors_correctable
     *          `t' defines the min hamming distance the resulting code must have
     *      * int n_data_bits
     *          determines the `m' and `n' parameters- the smallest length code of `t' 
     *          errors correctable is selected
     *      * int permutation
     *          specifies which primitive polynomial in the relevant GF(2^m) is used 
     *          (in the event there exists more than one). Permutation is 'modded' by the
     *          total number of primitive polynomials there are
     */
    class bch : public ecc_code
    {
        const static bool print_errors = false;
        bool initialized;
        int permutation;

        /*
         * Code parameters (based on the bch3.c implementation by Robert Morelos-Zaragoza)
         * 
         * m = order of the Galois field GF(2**m) 
         * n = 2**m - 1 = size of the multiplicative group of GF(2**m)
         * length = length of the BCH code codeword
         * k = n - deg(g(x)) = number of data bits theoretically storable in this code
         * t = error correcting capability (max. no. of errors the code corrects)
         * hd = 2*t + 1 = designed min. distance = no. of consecutive roots of g(x) + 1
         * alpha_to [] = log table of GF(2**m) 
         * index_of[] = antilog table of GF(2**m)
         * g[] = coefficients of the generator polynomial, g(x)
         * primitive_polynomial[] = coefficients of a primitive polynomial used to generate GF(2**m)
         */
        int m, n, length, k, t, hd;
        int n_data_bits; // the actual number of data bits that we make use of as requested at init time (i.e., potentially truncated)

        Eigen::Matrix< int, Eigen::Dynamic, 1 > alpha_to;
        Eigen::Matrix< int, Eigen::Dynamic, 1 > index_of;

        Eigen::Matrix< int, Eigen::Dynamic, 1 > g;
        Eigen::Matrix< int, Eigen::Dynamic, 1 > primitive_polynomial;
    public:
        std::string polynomial_to_str(Eigen::Matrix< int, Eigen::Dynamic, 1 > &p);

        bch(int permutation, int desired_n_data_bits, int desired_n_correctable_errs) : 
            initialized(false), t(desired_n_correctable_errs)
        {
            n_data_bits = desired_n_data_bits;
            t = desired_n_correctable_errs;
            hd = 2 * t + 1;

            // look for a BCH code that will work for these parameters- no such code is guaranteed to exist!
            if(find_valid_bch_params(permutation, desired_n_data_bits, desired_n_correctable_errs
                , k, length, m, primitive_polynomial, alpha_to, index_of, g))
            {
                printf("[ERROR] unable to determine a valid BCH code with k: %d, t: %d\n", desired_n_data_bits, desired_n_correctable_errs);
                return;
            }
            n = (1 << m) - 1; // size of the GF (primitive -> base == 2)

            if(g_verbosity > 0)
                printf("[INFO] Found usable BCH code of perm: %d, m: %d, n: %d k: %d t: %d (%d data + %d parity = %d code bits), g(x): %s\n"
                    , permutation, m, n, k, t, get_n_data_bits(), get_n_code_bits() - get_n_data_bits(), get_n_code_bits(), polynomial_to_str(g).c_str());

            // save internal state
            this->permutation = permutation;

            // only set on successful init
            initialized = true;
        }
        ~bch(void) {}

        bool ready(void) const
        {
            return initialized;
        }

        int correction_capability(void) const
        {
            return t;
        }

        int get_n_data_bits(void) const
        {
            return n_data_bits; // NOT k - k is the max number of data bits this code can support
        }

        int get_n_code_bits(void) const
        {
            return n - k + n_data_bits;
        }

        int get_permutation(void) const
        {
            return permutation;
        }

        std::string name(void) const
        {
            std::stringstream ss;
            ss << einsim::bch::static_name() << " (m: " << m << ", n: " << n << ", k: " << k << ", t: " << t << ") with ";
            ss << "#errors correctable: " << this->correction_capability() 
                << " (permutation: " << this->get_permutation()
                << ", n_data_bits: " << this->get_n_data_bits() 
                << ", n_code_bits: " << this->get_n_code_bits() << ")";
            return ss.str();            
        }

        std::string name_short(void) const
        {
            std::stringstream ss;
            ss << einsim::bch::static_name_short() << ": ";
            ss << "p:" << permutation << " t:" << t << " k:" << k << " n:" << n << " m:" << m;
            return ss.str();            
        }

        enum einsim::ecc_scheme get_scheme(void) const 
        {
            switch(this->correction_capability())
            {
                case 1: return einsim::ES_BCH_T1;
                case 2: return einsim::ES_BCH_T2;
                case 3: return einsim::ES_BCH_T3;
                default:
                    printf("[ERROR] unhandled correction capability for BCH code: %d\n", this->correction_capability());
                    assert(0 && "unhandled correction_capability for BCH code");
            }
        }

        static std::string static_name()
        {
            return std::string("BCH Code");
        }

        static std::string static_name_short()
        {
            return std::string("BCH");
        }

        int to_json(std::string &json) const { assert(0 && "unimplemented!"); }

        Eigen::Matrix< ET, Eigen::Dynamic, 1 > encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word) const;
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word) const;

        static void submit_tests(thread_pool &tp, enum test_mode mode);
    };
}

#endif /* BCH_CODE_H */
