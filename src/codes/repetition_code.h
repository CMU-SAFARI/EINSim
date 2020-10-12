/**
 * @file repetition_code.h
 *
 * @brief Encoding/decoding for n-repetition error correction codes
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef REPETITION_CODE_H
#define REPETITION_CODE_H

#include "Eigen/Eigen"
#include "supporting_routines.h"
#include "ecc_code.h"

namespace einsim
{
    /**
     * @brief implements an n-repetition code
     * 
     * The Repetition code is defined by the:
     *  * n_data_bits - the length of the dataword
     *  * n repetitions - effectively the correction capability of the code
     *  * permutation - a sort of ``random seed'' value that randomizes the
     *    locations of data bits within the codeword
     */
    class repetition : public ecc_code
    {
        bool initialized;
        int permutation;
        int n_data_bits;
        int n_reps;
        Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > bit_mapping_unshuffled;
        Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > bit_mapping;

    public:
        repetition(int permutation, int n_data_bits, int n_reps);
        ~repetition(void) {}

        bool ready(void) const
        {
            return initialized;
        }

        int correction_capability(void) const
        {
            return (n_reps - 1) / 2;
        }

        int get_n_data_bits(void) const
        {
            return n_data_bits;
        }

        int get_n_code_bits(void) const
        {
            return n_data_bits * n_reps;
        }

        int get_permutation(void) const
        {
            return permutation;
        }

        std::string name(void) const
        {
            std::stringstream ss;
            ss << einsim::repetition::static_name() << " with ";
            ss << "#errors correctable: " << this->correction_capability() << " (permutation: " << this->get_permutation() << ", n_reps: " << n_reps << ", n_data_bits: " << n_data_bits << ")";
            return ss.str();            
        }

        std::string name_short(void) const
        {
            std::stringstream ss;
            ss << einsim::repetition::static_name_short() << ": ";
            ss << "p:" << this->get_permutation() << " t:" << this->correction_capability() << " k:" << this->get_n_data_bits() << " n:" << this->get_n_code_bits();
            return ss.str();            
        }

        enum einsim::ecc_scheme get_scheme(void) const 
        {
            switch(this->correction_capability())
            {
                case 1: return einsim::ES_REPETITION_T1;
                case 2: return einsim::ES_REPETITION_T2;
                case 3: return einsim::ES_REPETITION_T3;
                default:
                    printf("[ERROR] unhandled correction capability for repetition code: %d\n", this->correction_capability());
                    assert(0 && "unhandled correction_capability for repetition code");
            }
        }

        static std::string static_name()
        {
            return std::string("Repetition Code");
        }

        static std::string static_name_short()
        {
            return std::string("REP");
        }

        int to_json(std::string &json) const { assert(0 && "unimplemented!"); }

        Eigen::Matrix< ET, Eigen::Dynamic, 1 > encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word) const;
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word) const;

        static void submit_tests(thread_pool &tp, enum test_mode mode);
    };
}

#endif /* REPETITION_CODE_H */