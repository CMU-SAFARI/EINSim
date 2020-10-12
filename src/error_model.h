/**
 * @file error_model.h
 *
 * @brief Data structures and routines used to represent different DRAM error models
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#ifndef ERROR_INJECTOR_H
#define ERROR_INJECTOR_H

#include <string>
#include <algorithm>
#include <random>

/* libraries */
#include "Eigen/Eigen"

/* project includes */
#include "word_generator.h"
#include "supporting_routines.h"

namespace einsim
{
    /* forward declaration */
    class error_model_descriptor;

    /** @brief enumerates different error distributions that the error injector can use */
    enum error_model
    {
          EM_NORMAL /**< normal cell that never fails */
        , EM_UNIFORM_RANDOM /**< uniform-random errors */
        , EM_DATA_RETENTION /**< uniform-random errors when a cell is programmed to the CHARGED state */
        , EM_DATA_RETENTION_NOISY /**< uniform-random errors when a cell is programmed to the CHARGED state with some random noise */
        , EM_STUCK_AT /**< cell that is always read as a fixed value no matter what is written */

        , EM_UNKNOWN
    };

    /* routines for converting between the string and enum representations of error distributions */
    enum einsim::error_model str_to_enum_error_model(const std::string &str); /**< @brief converts error_model string to enumeration */
    const std::string &enum_to_str_error_model(enum einsim::error_model em); /**< @brief converts error_model enumeration to string */
    std::string get_all_possible_error_models(void); /**< @brief returns a comma-separated string of all error_model values */
    
    /* routines to parse the fault injection mask from the CLI */
    // Eigen::Matrix< ET, Eigen::Dynamic, 1 > custom_fault_injection_mask_to_vector(const std::string &str); /**< @brief converts string representation of an error mask to an Eigen::Matrix */
    // std::string fault_injection_mask_to_str(const Eigen::Matrix< ET, Eigen::Dynamic, 1 > &str); /**< @brief converts an Eigen::Matrix representation of a fault_injection_mask to a string */
     
    /**
     * @brief probabilistically injects errors into a word with the given parameters
     * 
     * Injects errors according to the specified error model. This is
     * the primary error injection routine that simulates pre-correction error
     * characteristics.
     * 
     * @param word the word in which to inject errors
     * @param emd the error model 
     * @param dp the programmed data pattern to assume when injecting errors
     * @param tacs the true-/anti-cell state of the word that we're injecting errors in
     */
    void inject(
          Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
        , enum einsim::data_pattern dp
        , enum einsim::true_anti_cell_state tacs
        , const std::vector< einsim::error_model_descriptor * > &emd
    );
    
    /**
     * @brief injects a precise number of errors into a word with the given parameters
     * 
     * Generally used for testing and debugging, this function allows
     * injecting a precise number of errors into a word, though the actual
     * location of injected errors will vary according to the provided
     * parameters.
     * 
     * @param word the word in which to inject errors
     * @param em the error model 
     * @param cd the true-/anti-cell distribution to assume when injecting errors
     * @param dp the programmed data pattern to assume when injecting errors
     * @param n_errors the number of errors to inject
     */
    void inject_n(
          Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
        , enum error_model em
        , enum true_anti_cell_distribution cd
        , enum data_pattern dp
        , int n_errors
    );

    //
    // factory functions to produce the appropriate error model from a json file
    //
    /**
     * @brief extracts an error_model_descriptor object from a JSON configuration file
     * 
     * @param json_cfg_file JSON configuration filename
     * @return std::vector< std::vector< einsim::error_model_descriptor * > > list of extracted error model descriptors
     */
    std::vector< std::vector< einsim::error_model_descriptor * > > error_model_descriptors_from_json(const std::string &json_cfg_file);

    /**
     * @brief generates an error model descriptor object from a set of model parameters
     * 
     * @param em erorr model type to use
     * @param model_params list of model parametrs
     * @return einsim::error_model_descriptor* error model descriptor object
     */
    einsim::error_model_descriptor *error_model_descriptor_from_params(enum error_model em, const std::vector< std::string > &model_params);

    /**
     * @brief convert a list of error_model_descriptor objects to a printable string
     * 
     * @param emd_vec list of erorr model descriptors
     * @return std::string string representation of the error model descriptors
     */
    std::string error_model_descriptor_vec_to_str(const std::vector< einsim::error_model_descriptor * > &emd_vec);
    
    /**
     * @brief constructs the cartesian product of the set of provided error models
     * 
     * @param error_models cartesian product of the emds_per_bit
     * @param emds_per_bit per-bit error models to take the product over
     */
    void construct_cartesian_product_of_per_bit_error_models(
          std::vector< std::vector< einsim::error_model_descriptor * > > &error_models
        , const std::vector< std::vector< einsim::error_model_descriptor * > > &emds_per_bit);


    /**
     * @brief general base descriptor class for different types of error models that EINSim supports
     */
    class error_model_descriptor
    {
    protected:
        enum error_model em;

    public:
        error_model_descriptor(enum error_model em) : em(em) {}
        virtual ~error_model_descriptor() {};
        
        static int get_n_model_params(enum error_model em);
        
        enum error_model get_error_model(void) const { return em; }
        // std::string get_error_model_string(void) const { return enum_to_str_error_model(em); }
        
        virtual bool evaluate(bool data, bool is_true_cell) const = 0;
        virtual const std::string to_str(void) const = 0;
    };
    

    /**
     * @brief representation of a bit that never experiences error
     */
    class error_model_descriptor_normal : public error_model_descriptor
    {
    public:
        /**
         * @brief state parameters necessary to keep track of... well... nothing :)
         */
        struct params
        {

            std::string to_str(void) const
            {
                return "";
            }

            params() {}
            ~params() {}
        } p;

    protected:
        
    public:
        error_model_descriptor_normal() : 
              error_model_descriptor(einsim::EM_NORMAL)
        {
        }
        
        static int get_n_model_params(void)
        {
            return 0;
        }

        bool evaluate(bool data, bool is_true_cell) const
        {
            return data;
        }

        const std::string to_str(void) const
        {
            return enum_to_str_error_model(this->get_error_model()) + "(" + p.to_str() + ")";
        }
    };

    /**
     * @brief models uniform-random errors irrespective of the cell type
     */
    class error_model_descriptor_uniform_random : public error_model_descriptor
    {
    public:
        /**
         * @brief state parameters necessary to keep track of the uniform-random error distribution
         */
        struct params
        {
            float p;

            std::string to_str(void) const
            {
                return std::string("p:") + std::to_string(p);
            }

            params(float error_rate) : p(error_rate)  {}
            ~params() {}
        } p;

    protected:
        mutable std::default_random_engine rng;
        mutable std::bernoulli_distribution dist;
        
    public:
        error_model_descriptor_uniform_random(float error_rate) : 
              error_model_descriptor(einsim::EM_UNIFORM_RANDOM)
            , p(error_rate)
            , dist(error_rate)
        {
            rng.seed(std::random_device{}());
        }
        
        static int get_n_model_params(void)
        {
            return 1;
        }

        bool evaluate(bool data, bool is_true_cell) const
        {
            return dist(rng) ? !data : data;
        }

        const std::string to_str(void) const
        {
            return enum_to_str_error_model(this->get_error_model()) + "(" + p.to_str() + ")";
        }
    };

    /**
     * @brief uniform-random probability from the CHARGED to the DISCHARGED state
     */
    class error_model_descriptor_data_retention : public error_model_descriptor
    {
    public:
        /**
         * @brief state parameters necessary to keep track of the uniform-random error distribution
         */
        struct params
        {
            float p;

            std::string to_str(void) const
            {
                return std::string("p:") + std::to_string(p);
            }

            params(float error_rate) : p(error_rate)  {}
            ~params() {}
        } p;


    protected:
        mutable std::default_random_engine rng;
        mutable std::bernoulli_distribution dist;
        

    public:
        error_model_descriptor_data_retention(float error_rate) : 
              error_model_descriptor(einsim::EM_DATA_RETENTION)
            , p(error_rate)
            , dist(error_rate)
        {
            rng.seed(std::random_device{}());            
        }
        
        static int get_n_model_params(void)
        {
            return 1;
        }

        bool evaluate(bool data, bool is_true_cell) const
        {
            bool can_fail = data == is_true_cell;
            if(can_fail)
                return dist(rng) ? !data : data;
            return data;
        }

        const std::string to_str(void) const
        {
            return enum_to_str_error_model(this->get_error_model()) + "(" + p.to_str() + ")";
        }
    };

    /**
     * @brief uniform-random probability from the CHARGED to the DISCHARGED state with a noise component
     */
    class error_model_descriptor_data_retention_noisy : public error_model_descriptor
    {
    public:
        /**
         * @brief state parameters necessary to keep track of the uniform-random error distribution with additional noise
         */
        struct params
        {
            float p; // error probability
            float n; // noise ratio

            std::string to_str(void) const
            {
                return std::string("p:") + std::to_string(p) + std::string(" n:") + std::to_string(n);
            }

            params(float error_rate, float noise_ratio) : p(error_rate), n(noise_ratio) {}
            ~params() {}
        } p;


    protected:
        mutable std::default_random_engine rng;
        mutable std::bernoulli_distribution error_dist;
        mutable std::bernoulli_distribution noise_dist;

    public:
        error_model_descriptor_data_retention_noisy(float error_rate, float noise_ratio) : 
              error_model_descriptor(einsim::EM_DATA_RETENTION_NOISY)
            , p(error_rate, noise_ratio)
            , error_dist(error_rate)
            , noise_dist(noise_ratio)
        {
            rng.seed(std::random_device{}());            
        }
        
        static int get_n_model_params(void)
        {
            return 2;
        }

        bool evaluate(bool data, bool is_true_cell) const
        {
            bool can_fail = data == is_true_cell;
            if(can_fail)
                data = error_dist(rng) ? !data : data;
            return noise_dist(rng) ? !data : data;
        }

        const std::string to_str(void) const
        {
            return enum_to_str_error_model(this->get_error_model()) + "(" + p.to_str() + ")";
        }
    };

    /**
     * @brief cell that is fixed at a given value no matter what is written
     */
    class error_model_descriptor_stuck_at : public error_model_descriptor
    {
    public:
        /**
         * @brief state parameters necessary to keep track of the stuck-at error model
         */
        struct params
        {
            bool v;

            std::string to_str(void) const
            {
                return std::string("v:") + std::to_string(v);
            }

            params(bool stuck_at_val) : v(stuck_at_val)  {}
            ~params() {}
        } p;

    protected:
        
    public:
        error_model_descriptor_stuck_at(bool stuck_at_val) : 
              error_model_descriptor(einsim::EM_STUCK_AT)
            , p(stuck_at_val)
        {
        }
        
        static int get_n_model_params(void)
        {
            return 1;
        }

        bool evaluate(bool data, bool is_true_cell) const
        {
            return p.v;
        }

        const std::string to_str(void) const
        {
            return enum_to_str_error_model(this->get_error_model()) + "(" + p.to_str() + ")";
        }
    };
}

#endif /* ERROR_INJECTOR_H */
