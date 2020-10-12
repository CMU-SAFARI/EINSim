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

/* library routines */
#include "Eigen/Eigen"
#include "rapidjson/document.h"

/* project includes */
#include "supporting_routines.h"
#include "ecc_code.h"

namespace einsim
{
    /**
     * @brief implements a single-error correction (SEC) Hamming code
     * 
     * The Hamming code is defined by the:
     *      * n_data_bits - the length of the dataword
     *      * permutation - a sort of ``random seed'' value that randomizes
     *        the location of data bits in the codeword
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
        std::vector< int > data_bit_indices;
        std::vector< int > parity_bit_indices;

    public:
        hamming(int permutation, int n_data_bits);
        hamming(const rapidjson::Document &d, const std::string &cfg_file_name); // initialize from configuration file
        ~hamming(void) {}

        bool ready(void) const
        {
            return initialized;
        }

        int correction_capability(void) const
        {
            return 1;
        }

        int get_n_data_bits(void) const
        {
            return nd;
        }

        int get_n_code_bits(void) const
        {
            return nd + np;
        }

        int get_permutation(void) const
        {
            return permutation;
        }

        std::string name(void) const
        {
            std::stringstream ss;
            ss << einsim::hamming::static_name() << " with ";
            ss << "#errors correctable: " << this->correction_capability() 
                << " (permutation: " << this->get_permutation() 
                << ", n_data_bits: " << nd 
                << ", n_parity_bits: " << np << ")";
            return ss.str();            
        }

        std::string name_short(void) const
        {
            std::stringstream ss;
            ss << einsim::hamming::static_name_short() << ": ";
            ss << "p:" << this->get_permutation() << " t:" << this->correction_capability() << " k:" << this->get_n_data_bits() << " n:" << this->get_n_code_bits();
            return ss.str();            
        }

        enum einsim::ecc_scheme get_scheme(void) const 
        {
            return einsim::ES_HAMMING_SEC;
        }

        static std::string static_name()
        {
            return std::string("Hamming SEC Code");
        }

        static std::string static_name_short()
        {
            return std::string("HSC");
        }

        uint64_t compute_uid(void) const;
        int to_json(std::string &json) const;

        Eigen::Matrix< ET, Eigen::Dynamic, 1 > encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word) const;
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word) const;

        static void submit_tests(thread_pool &tp, enum test_mode mode);
    };
}

#endif /* HAMMING_CODE_H */