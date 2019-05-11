/**
 * @file hamming_code.cpp
 *
 * @brief Implementation of the Hamming SEC code
 * 
 * Note: these functions must be threadsafe since they'll be called from different threads.
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
/* stdlib */
#include <set>
#include <string>
#include <vector>
#include <iostream>

/* libraries */
#include "Eigen/Eigen"
#include "libtp.h"

/* project includes */
#include "supporting_routines.h"
#include "ecc_code.h"
#include "hamming_code.h"

einsim::hamming::hamming(int permutation, int n_data_bits)
	: initialized(false), permutation(permutation), nd(n_data_bits) 
{
	if(n_data_bits <= 0)
	{
		printf_both("ERROR: Invalid number of data bits: %d\n", n_data_bits);
		exit(-1);
	}

	// compute n_parity_bits
	np = 0;
	while((1 << np) < (np + nd + 1))
		np++;
	assert(np == (int) std::ceil(std::log2(nd + np + 1)));

	// std::cout << "Creating hamming code of " << nd << " data bits, " << np << " parity bits" << std::endl;

	// generate a set of each non-Po2 value
	std::vector< int > syn_values;
	for(int i = 0; i < (1 << np); i++)
		if(i & (i - 1))
			syn_values.push_back(i);

	// select a random subset of (nd) items
	std::seed_seq seed{permutation};
	std::mt19937 rng(seed);
	std::vector< int > indices(syn_values.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::shuffle(indices.data(), indices.data() + indices.size(), rng);

	// construct the parity-check matrix in form [P^t | I_{n-k}]
    parity_check = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(np, nd + np);
	assert(parity_check.rows() == np);
	assert(parity_check.cols() == np + nd);			
	for(int col = 0; col < nd; col++) // inject syndromes
		for(int row = 0; row < np; row++)
			parity_check(row, col) = (syn_values[indices[col]] >> row) & 1;	
	parity_check.block(0, nd, np, np) 
		= Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic >::Identity(np, np);

	// construct the generator matrix
	generator = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(nd + np, nd);
	assert(generator.rows() == nd + np);
	assert(generator.cols() == nd);
	for(int row = 0; row < np; row++)
			generator.row(nd + row) = parity_check.row(row).segment(0, nd);
	generator.block(0, 0, nd, nd) 
		= Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic >::Identity(nd, nd);
	
	// construct the degenerator matrix (rows represent order of bits to choose)
	degenerator = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(nd, nd + np);
	degenerator.block(0, 0, nd, nd) 
		= Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic >::Identity(nd, nd);

	// std::cout << parity_check << std::endl;
	// std::cout << generator << std::endl;
	
	// shuffle the cols of H, rows of G
	// assert(parity_check.rows() == np);
	// assert(parity_check.cols() == np + nd);

	// randomize the rows according to the 'permutation'
	Eigen::PermutationMatrix< Eigen::Dynamic, Eigen::Dynamic > perm(nd + np);
	perm.setIdentity();
	std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(), rng);

	// permute columns of matrices
	// A_perm = A * perm; // permute columns


	initialized = true;
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::hamming::encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word)
{
	assert(data_word.size() == nd);
	return MOD2_EIGEN_VEC(generator * data_word);
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::hamming::decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word)
{
	assert(code_word.size() == nd + np);
	Eigen::Matrix< ET, Eigen::Dynamic, 1 > syndrome = MOD2_EIGEN_VEC(parity_check * code_word);

	/* find syndrome in matrix H */
	for(int col = 0; col < nd + np; col++)
	{
		if(parity_check.col(col).isApprox(syndrome))
		{
			code_word[col] ^= 1;
			break;
		}
	}
	
	return degenerator * code_word;
}

// tests whatever configurations we need to test- this is ECC type-specific so it goes here
void einsim::hamming::submit_tests(thread_pool &tp, enum test_mode mode)
{
	fprintf(g_output_file, "Testing %s\n", einsim::hamming::static_name().c_str());

	// generic lambda, operator() is a template with one parameter
	auto lambdafunc = [](int tid, int niter, int perm, int n_db, thread_pool *tp) 
		{
	        einsim::ecc_code *ec = new einsim::hamming(perm, n_db);
	        for(int iterations = 0; iterations < niter; iterations++)
	            tp->add(einsim::test_thread, 0 /* priority */, ec);
    	};

    if(mode == einsim::TM_SLOW)
    {
	    for(int perm = 0; perm < 10; perm++)
	    	for(int n_db = 1; n_db < 1000; n_db <<= 1)
	    	{
	    		if(n_db > 1)
	    			tp.add(lambdafunc, 1, 100, perm, n_db - 1, &tp);
	    		tp.add(lambdafunc, 1, 100, perm, n_db, &tp);
	    		tp.add(lambdafunc, 1, 100, perm, n_db + 1, &tp);
	    	}
    }
	else if(mode == einsim::TM_FAST)
	{
	    for(int perm = 0; perm < 10; perm++)
		    for(int n_db : std::set< int >({1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 63, 64, 65, 127, 128, 129, 255, 256}))
	    		tp.add(lambdafunc, 1, 1, perm, n_db, &tp);
	}
	else
	{
		printf_both("Invalid test mode: %d\n", mode);
		exit(-1);
	}
}