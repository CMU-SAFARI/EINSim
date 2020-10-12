/**
 * @file repetition_code.h
 *
 * @brief Implementation of n-repetition code
 * 
 * Note: these functions must be threadsafe since they'll be called from different threads.
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
/* stdlib */
#include <iostream>
#include <set>
#include <string>
#include <vector>

/* libraries */
#include "libtp.h"
#include "Eigen/Eigen"

/* project */
#include "supporting_routines.h"
#include "ecc_code.h"
#include "repetition_code.h"

einsim::repetition::repetition(int permutation, int n_data_bits, int n_reps) 
    : initialized(false), permutation(permutation), n_data_bits(n_data_bits), n_reps(n_reps) 
{
    if(n_data_bits <= 0)
    {
        printf("ERROR: Invalid number of data bits: %d\n", n_data_bits);
        exit(-1);
    }
    if((n_reps & 1) == 0)
    {
        printf("ERROR: Invlid number of repetitions: %d: ambiguous decoding for even repetitions\n", n_reps);
        exit(-1);
    }

    // std::cout << "Creating repetition code of " << n_data_bits << " data bits, " << n_reps << " repetition" << std::endl;

    // Encoding: Ax = b -> [16*Nx16][16] = [16*N]
    // Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > bit_mapping_unshuffled;
    bit_mapping_unshuffled.resize(n_data_bits * n_reps /* rows */, n_data_bits /* cols */);
    for(int n = 0; n < n_data_bits; n++)
    {
        for(int r = 0; r < n_reps; r++)
        {
            Eigen::Matrix< ET, 1, Eigen::Dynamic > row_vector(n_data_bits);
            for(int c = 0; c < n_data_bits; c++)
                row_vector[c] = !!(c == n);
            bit_mapping_unshuffled.row(n * n_reps + r) = row_vector;
        }
    }
    // std::cout << "Mapping: " << std::endl << bit_mapping << std::endl;

    // randomize the rows according to the 'permutation'
    Eigen::PermutationMatrix< Eigen::Dynamic, Eigen::Dynamic > perm(n_data_bits * n_reps);
    perm.setIdentity();

    std::seed_seq seed{permutation};
    std::mt19937 rng(seed);
    std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(), rng);

    Eigen::Matrix< ET, Eigen::Dynamic, 1 > data(n_data_bits);
    for(int n = 0; n < n_data_bits; n++)
        data[n] = rand() & 1;
    // std::cout << "Data: " << data << std::endl;
    // std::cout << "Codeword: " << bit_mapping * data << std::endl;
    // Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_encoded = bit_mapping * data;
    bit_mapping = perm * bit_mapping_unshuffled;
    // Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_encoded_shuffled = bit_mapping * data;
    // assert(data_encoded.sum() == data_encoded_shuffled.sum());


    // std::cout << "Codeword shuffled: " << bit_mapping * data << std::endl;
    // std::cout << "Mapping shuffled: " << std::endl << bit_mapping << std::endl;

    initialized = true;
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::repetition::encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word) const
{
    // Encoding: Ax = b -> [16*Nx16][16] = [16*N]
    return bit_mapping * data_word;
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::repetition::decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word) const
{
    // Decoding: x = A^-1 b
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word = bit_mapping.transpose() * code_word;
    return (data_word * 2) / (n_reps + 1);
}

// tests whatever configurations we need to test- this is ECC type-specific so it goes here
void einsim::repetition::submit_tests(thread_pool &tp, enum test_mode mode)
{
    std::cout << "Testing " << einsim::repetition::static_name() << std::endl;

    // generic lambda, operator() is a template with one parameter
    auto lambdafunc = [](int tid, int niter, int ep, int n_db, int nr, thread_pool *tp) 
        {
            einsim::ecc_code *ec = new einsim::repetition(ep, n_db, nr);
            for(int iterations = 0; iterations < niter; iterations++)
                tp->add(einsim::test_thread, 0 /* priority */, ec);
        };

    if(mode == einsim::TM_SLOW)
    {
        for(int ep = 0; ep < 10; ep++)
            for(int n_db : std::set< int >({1, 2, 3, 4, 7, 8, 15, 16, 31, 32
                , 63, 64, 127, 128, 255, 256, 511, 512}))
                for(int nr = 3; nr <= 11; nr += 2)
                    tp.add(lambdafunc, 1, 100, ep, n_db, nr, &tp);
    }
    else if(mode == einsim::TM_FAST)
    {
        for(int ep = 0; ep < 2; ep++)
            for(int n_db : std::set< int >({1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 64, 128, 256}))
                for(int nr = 3; nr <= 9; nr += 2)
                    tp.add(lambdafunc, 1, 1, ep, n_db, nr, &tp);
    }
    else
        assert(0 && "Invalid test mode");
}