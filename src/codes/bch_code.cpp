/**
 * @file bch_code.cpp
 *
 * @brief Implementation of the BCH code
 *     
 * This implementation relies heavily on the helper functions in bch_helpers.*
 * to do the heavy lifting.
 * 
 * Note: these functions must be threadsafe since they'll be called from different threads.
 *  
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
/* stdlib */
#include <string>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>

#include "supporting_routines.h"
#include "ecc_code.h"
#include "codes/bch_helpers.h"
#include "codes/bch_code.h"

std::string einsim::bch::polynomial_to_str(Eigen::Matrix< int, Eigen::Dynamic, 1 > &p)
{
    std::stringstream ss;
    ss << "0b";
    for(int i = p.size() - 1; i >= 0; i--) 
        ss << (p[i] & 1);
    ss << ", 0o";
    for(int digit = (p.size() + 2) / 3 - 1; digit >= 0; digit--) 
    {
        uint8_t oct = 0;
        for(int i = 2; i >= 0; i--)
        {
            int pidx = digit * 3 + i;
            if(pidx < p.size())
                oct |= p[pidx] << i;
        }
        ss << (oct & 0x07);
    }
    ss << ", 0x";
    for(int digit = (p.size() + 3) / 4 - 1; digit >= 0; digit--) 
    {
        uint8_t oct = 0;
        for(int i = 2; i >= 0; i--)
        {
            int pidx = digit * 4 + i;
            if(pidx < p.size())
                oct |= p[pidx] << i;
        }
        ss << (oct & 0xf);
    }
    return ss.str();
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::bch::encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word)
{
    assert(data_word.size() == get_n_data_bits());

    // 0-pad up to k for symbol computation
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_padded(k);
    if(k - n_data_bits)
        data_word_padded << data_word, Eigen::Matrix< ET, Eigen::Dynamic, 1 >::Zero(k - n_data_bits);
    else
        data_word_padded = data_word;
    assert(data_word_padded.size() == k);

    // call the Robert Morelos-Zaragoza's decoding routine, which is modified
    // to try and fix as many errors as it can
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > parity_bits = bch_encode(length, k, data_word_padded, g);

    Eigen::Matrix< ET, Eigen::Dynamic, 1 > codeword(get_n_code_bits());
    codeword << parity_bits, data_word; // this is a systematic code
    assert(codeword.size() == get_n_code_bits());
    return codeword;
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::bch::decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word)
{
    assert(code_word.size() == get_n_code_bits());

    // 0-pad up to k for decoding
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word_padded(n);
    if(k - n_data_bits)
        code_word_padded << code_word, Eigen::Matrix< ET, Eigen::Dynamic, 1 >::Zero(k - n_data_bits);
    else
        code_word_padded = code_word;
    assert(code_word_padded.size() == n);

    // call the Robert Morelos-Zaragoza's decoding routine, which is modified
    // to try and fix as many errors as it can
    bch_decode(length, t, n, code_word_padded, alpha_to, index_of);

    // unpad the word to extract the data bits alone
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_unpadded = code_word_padded.bottomRows(k).topRows(n_data_bits);
    assert(data_word_unpadded.size() == get_n_data_bits());
    return data_word_unpadded;
}

// tests whatever configurations we need to test- this is ECC type-specific so it goes here
void einsim::bch::submit_tests(thread_pool &tp, enum test_mode mode)
{
    std::cout << "Testing " << einsim::bch::static_name() << std::endl;

    // generic lambda, operator() is a template with one parameter
    auto lambdafunc = [](int tid, int niter, int ecc_perm, int n_desired_data_bits, int nerrs_correctable, thread_pool *tp) 
        {
            einsim::ecc_code *ec = new einsim::bch(ecc_perm, n_desired_data_bits, nerrs_correctable);
            for(int iterations = 0; iterations < niter; iterations++)
                tp->add(einsim::test_thread, 0 /* priority */, ec);
        };

    if(mode == einsim::TM_SLOW)
    {
        // test a sweep over the configuration parameters
        for(int ep = 0; ep < 10; ep++)
            for(int code_len : std::set< int >({1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 64, 128, 256}))
                for(int nerrs_correctable = 3; nerrs_correctable <= 9; nerrs_correctable += 2)
                    tp.add(lambdafunc, 1, 100, ep, code_len, nerrs_correctable, &tp);
    }
    else if(mode == einsim::TM_FAST)
    {
        // test specific tuples
        std::vector< std::tuple< int, int, int > > test_tuples = 
        {
              std::make_tuple(0, 128, 1)
            , std::make_tuple(0, 128, 2)
            , std::make_tuple(0, 128, 3)
            , std::make_tuple(0, 128, 4)
            , std::make_tuple(0, 128, 5)
            , std::make_tuple(0, 128, 6)
            , std::make_tuple(0, 128, 7)
        };
        for(const auto &tuple : test_tuples)
            tp.add(lambdafunc
              , 1 /* priority */
              , 100 /* niter */
              , std::get< 0 >(tuple) /* ecc_perm */
              , std::get< 1 >(tuple) /* n_desired_data_bits */
              , std::get< 2 >(tuple) /* nerrs correctable */
              , &tp);
    }
    else
        assert(0 && "Invalid test mode");
}