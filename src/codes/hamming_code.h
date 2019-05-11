/**
 * @file hamming_code.h
 *
 * @brief Encoding/decoding for Hamming single-error correction codes
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef HAMMING_CODE_H
#define HAMMING_CODE_H

#include "Eigen/Eigen"
#include "supporting_routines.h"
#include "ecc_code.h"

namespace einsim
{
	/**
	 * @brief implements a single-error correction (SEC) Hamming code
	 * 
	 * The Hamming code is defined by the:
	 * 		* n_data_bits - the length of the dataword
	 * 		* permutation - a sort of ``random seed'' value that randomizes
	 * 		  the location of data bits in the codeword
	 */
	class hamming : public ecc_code
	{
		bool initialized;
		int permutation;
		int nd;
		int np;
		Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > generator;
		Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > degenerator;
		Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > parity_check;

	public:
		hamming(int permutation, int n_data_bits);
		~hamming(void) {}

		bool ready(void)
		{
			return initialized;
		}

		int correction_capability(void)
		{
			return 1;
		}

		int get_n_data_bits(void)
		{
			return nd;
		}

		int get_n_code_bits(void)
		{
			return nd + np;
		}

		int get_permutation(void)
		{
			return permutation;
		}

		std::string name(void)
		{
			std::stringstream ss;
			ss << einsim::hamming::static_name() << " with ";
			ss << "#errors correctable: " << this->correction_capability() 
				<< " (permutation: " << this->get_permutation() 
				<< ", n_data_bits: " << nd 
				<< ", n_parity_bits: " << np << ")";
			return ss.str();			
		}

		std::string name_short(void)
		{
			std::stringstream ss;
			ss << einsim::hamming::static_name_short() << ": ";
			ss << "p:" << this->get_permutation() << " t:" << this->correction_capability() << " k:" << this->get_n_data_bits() << " n:" << this->get_n_code_bits();
			return ss.str();			
		}

		static std::string static_name()
		{
			return std::string("Hamming SEC Code");
		}

		static std::string static_name_short()
		{
			return std::string("HSC");
		}

		Eigen::Matrix< ET, Eigen::Dynamic, 1 > encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word);
		Eigen::Matrix< ET, Eigen::Dynamic, 1 > decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word);

		static void submit_tests(thread_pool &tp, enum test_mode mode);
	};
}

#endif /* HAMMING_CODE_H */