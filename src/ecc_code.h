/**
 * @file ecc_code.h
 *
 * @brief Definitions of the chief data structures used throughout EINSim
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef ECC_CODE_H
#define ECC_CODE_H

/* stdlib */
#include <string>

/* libraries */
#include "Eigen/Eigen"
#include "libtp.h"

/* project includes */
#include "supporting_routines.h"

namespace einsim
{	
	/* forward declarations */
	class ecc_code;
	void test_thread(int tid, einsim::ecc_code *code);

	/** @brief enumeration representing different test modes available for ECC schemes */
	enum test_mode
	{
		  TM_FAST /**< fast test mode that covers basic functionality */
		, TM_SLOW /**< detailed test mode with higher coverage of the ECC scheme's parameter space */

		, TM_UNKNOWN
	};

	/*  convenience functions for converting between the enumeration and string representations of test modes */
	std::string enum_to_str_test_mode(enum einsim::test_mode es); /**< @brief converts test mode enumeration to string */
	enum einsim::test_mode str_to_enum_test_mode(std::string str); /**< @brief converts test mode string to enumeration */
	std::string get_all_possible_test_modes(void); /**< @brief returns all valid test modes as a comma-separated string */

	/** @brief enumeration representing different types of ECC codes that are implemented in EINSim */
	enum ecc_scheme
	{
		  ES_REPETITION_T1
		, ES_REPETITION_T2
		, ES_REPETITION_T3
		, ES_HAMMING_SEC
		, ES_BCH_T1
		, ES_BCH_T2
		, ES_BCH_T3

		, ES_UNKNOWN
	};

	/* convenience functions for converting between the enumeration and string representation of ecc schemes */
	std::string enum_to_str_ecc_scheme(enum einsim::ecc_scheme es); /**< @brief converts ecc scheme enumeration to string */
	std::string get_all_possible_ecc_schemes(void); /**< @brief returns all valid ecc schemes as a comma-separated string */

	/**
	 * @brief virtual base class representing the interface for an ECC code implementation 
	 */
	class ecc_code
	{
	private:

	public:
		ecc_code() {};
		~ecc_code() {};

		virtual std::string name(void) = 0;
		virtual std::string name_short(void) = 0;

		virtual int correction_capability(void) = 0;
		virtual int get_n_data_bits(void) = 0;
		virtual int get_n_code_bits(void) = 0;
		virtual int get_permutation(void) = 0;
		virtual bool ready(void) = 0;

		virtual Eigen::Matrix< ET, Eigen::Dynamic, 1 > encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word) = 0;	
		virtual Eigen::Matrix< ET, Eigen::Dynamic, 1 > decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word) = 0;	
	};
	
	/**
	 * @brief a single testing worker thread that tests one ECC scheme
	 * 
	 * @param tid worker thread ID
	 * @param ec pointer to the ECC scheme to test
	 */
	void test_thread(int tid, einsim::ecc_code *ec);

	/**
	 * @brief parent testing routine that tests ECC schemes with the given test mode and worker testing function
	 * 
	 * @param test_func the worker thread function to invoke
	 * @param mode the test mode to use
	 * @param n_threads number of threads to parallelize over
	 */
	void test_ecc(void (*test_func)(thread_pool &tp, enum test_mode mode), enum test_mode mode, int n_threads);
}

#endif /* ECC_CODE_H */