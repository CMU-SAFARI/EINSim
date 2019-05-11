/**
 * @file debug.h
 *
 * @brief Declares routines that are useful for debugging ECC implementations
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef DEBUG_H
#define DEBUG_H

/**
 * @brief Minimum unit of parallelism for the debug_example
 * 
 * Simulates a parameterizable number of words for a fixed ECC code. Used as the 
 * minimum unit of parallelism for the debug_example code.
 * 
 * @param tid worker thread ID
 * @param ec pointer to the ECC code to simulate
 * @param n_words_to_simulate total number of words to simulate
 * @param dp the data pattern to simulate
 */
void debug_example_worker(int tid, einsim::ecc_code *ec, int n_words_to_simulate, enum einsim::data_pattern dp);

/**
 * @brief Example of a targeted testing/debugging routine
 * 
 * Provides a working example of a simpler, more controllable error-injection
 * loop that may be used for testing/debugging a new feature or ECC scheme. It
 * generally allows for a well-defined environment that can target specific
 * parameters as required.
 *  
 * @param n_worker_threads number of worker threads to use
 * @param n_words_to_simulate total number of words to simulate
 */
void debug_example(int n_worker_threads, int n_words_to_simulate);

#endif /* DEBUG_H */