/**
 * @file error_model.cpp
 *
 * @brief Definitions of routines for representing different per-bit error models
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>

/* libraries */
#include "Eigen/Eigen"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"

/* project includes */
#include "supporting_routines.h"
#include "error_model.h"

/* random distribution parameters used by the error injector */
static std::default_random_engine rng;
static std::bernoulli_distribution ei_bernoulli_0_5(0.5);
static bool error_injector_initialized = false;

static std::random_device rd_inject_errs;
static std::mt19937 rng_inject_errs(rd_inject_errs());

enum einsim::error_model einsim::str_to_enum_error_model(const std::string &str)
{
    static std::map< std::string, enum error_model > str_to_enum_error_model_map = 
    {
          {"NORMAL", EM_NORMAL}
        , {"UNIFORM_RANDOM", EM_UNIFORM_RANDOM}
        , {"DATA_RETENTION", EM_DATA_RETENTION}
        , {"DATA_RETENTION_NOISY", EM_DATA_RETENTION_NOISY}
        , {"STUCK_AT", EM_STUCK_AT}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_error_model_map.find(str_uppercase);
    if(search == str_to_enum_error_model_map.end())
        return EM_UNKNOWN;
    else
        return search->second;
}

const std::string &einsim::enum_to_str_error_model(enum einsim::error_model em)
{
    static std::map< enum error_model, std::string > enum_to_str_error_model_map = 
    {
          {EM_NORMAL, "NORMAL"}
        , {EM_UNIFORM_RANDOM, "UNIFORM_RANDOM"}
        , {EM_DATA_RETENTION, "DATA_RETENTION"}
        , {EM_DATA_RETENTION_NOISY, "DATA_RETENTION_NOISY"}
        , {EM_STUCK_AT, "STUCK_AT"}
        
        , {EM_UNKNOWN, "UNKNOWN"}
    };
    auto search = enum_to_str_error_model_map.find(em);
    if(search == enum_to_str_error_model_map.end())
        return enum_to_str_error_model_map.at(einsim::EM_UNKNOWN);
    else
        return search->second;
}

std::string einsim::get_all_possible_error_models(void)
{
    std::string ret;
    for(int s = einsim::EM_NORMAL; s < einsim::EM_UNKNOWN; s++)
        ret += einsim::enum_to_str_error_model((einsim::error_model)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}

void einsim::inject(
      Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
    , enum einsim::data_pattern dp
    , enum einsim::true_anti_cell_state tacs
    , const std::vector< einsim::error_model_descriptor * > &emd
)
{    
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > word_pre_injection(word);

    // initialize the random distribution parameters
    if(!error_injector_initialized)
    {
        rng.seed(std::random_device{}());
        error_injector_initialized = true;
    }

    // inject errors
    for(int i = 0; i < word.size(); i++)
    {
        bool data = word(i);
        bool is_true_cell;
        switch(tacs)
        {
            // all cells in this word are TRUE cells
            case einsim::TACS_ALL_TRUE:
                is_true_cell = true;
                break;

            // all cells in this word are ANTI-cells
            case einsim::TACS_ALL_ANTI:
                is_true_cell = false;
                break;

            // cells alternate T/A starting with T
            case einsim::TACS_ALT_T:
                is_true_cell = (i % 2) == 0;
                break;

            // cells alternate T/A starting with A
            case einsim::TACS_ALT_A:
                is_true_cell = (i % 2) == 1;
                break;

            default:
                printf_both("[ERROR] Unsupported true-/anti-cell state: %d\n", (int)tacs);
                assert(0 && "unsupported true-/anti-cell state");
        }
        
        if(emd.size() == 1)
            word(i) = emd.at(0)->evaluate(data, is_true_cell);
        else
            word(i) = emd.at(i)->evaluate(data, is_true_cell);
    }
}

// used only for testing/debugging, so we only handle the CHARGED and DISCHARGED cases
// and any cases that result in the same behavior (e.g., ALL_TRUE + ALL_ONES)
void einsim::inject_n(
      Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
    , enum error_model em
    , enum true_anti_cell_distribution cd
    , enum data_pattern dp
    , int n_errors
)
{
    if(!error_injector_initialized)
    {
        rng.seed(std::random_device{}());
        error_injector_initialized = true;
    }

    // Ideally, we need to create a mask of the bits that can actually fail
    // depending on the cd + dp. However, that's VERY slow. Instead, we only handle the cases
    // in which we effectively end up CHARGED or DISCHARGED
    bool is_charged = false;
    switch(cd)
    {
        case CD_ALL_TRUE_OR_ALL_ANTI:
        {
            switch(dp)
            {
                case DP_RANDOM:
                    assert(0 && "there's no way to accurately compute this - the given parameters may or may not be impossible. use something deterministic or implement this under your own assumptions");
                    break;

                case DP_ALL_ONES:
                {
                    bool is_true_cell = ei_bernoulli_0_5(rng) ? 0 : 1;
                    if(is_true_cell)
                        is_charged = true;
                    else
                        is_charged = false;
                    break;
                }

                case DP_CHARGED:
                    is_charged = true;
                    break;

                default:
                    printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
                    assert(0);
            }
            break;
        }

        case CD_ALL_TRUE:
        {
            switch(dp)
            {
                case DP_RANDOM:
                    assert(0 && "there's no way to accurately compute this - the given parameters may or may not be impossible. use something deterministic or implement this under your own assumptions");
                    break;

                case DP_ALL_ONES:
                    is_charged = true;
                    break;

                case DP_CHARGED:
                    is_charged = true;
                    break;

                default:
                    printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
                    assert(0);
            }
            break;
        }

        case CD_ALL_ANTI:
        {
            switch(dp)
            {
                case DP_RANDOM:
                    assert(0 && "there's no way to accurately compute this - the given parameters may or may not be impossible. use something deterministic or implement this under your own assumptions");
                    break;

                case DP_ALL_ONES:
                    is_charged = false;
                    break;

                case DP_CHARGED:
                    is_charged = true;
                    break;

                default:
                    printf_both("[ERROR] Unsupported data pattern: %d\n", (int)dp);
                    assert(0);
            }
            break;
        }

        case CD_COLSTRIPE_T:
        {
            switch(dp)
            {
                default:
                    printf_both("[ERROR] Unimplemented data pattern: %d\n", (int)dp);
                    assert(0);
            }
            break;
        }

        case CD_COLSTRIPE_A:
        {
            switch(dp)
            {
                default:
                    printf_both("[ERROR] Unimplemented data pattern: %d\n", (int)dp);
                    assert(0);
            }
            break;
        }

        default:
            printf_both("[ERROR] Unsupported true-/anti-cell distribution: %d\n", (int)cd);
            assert(0);
            break;
    }

    // check if it's even possible to satisfy the request
    if(is_charged)
    {
        assert(n_errors <= word.size() && "impossible to meet requests for the given configuration");
        assert(em == EM_UNIFORM_RANDOM && "unsupported error model requested!");

        // inject N errors uniform-randomly
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > mask(word.size());
        int i;
        for(i = 0; i < n_errors; i++)
            mask[i] = 1;
        for(; i < word.size(); i++)
            mask[i] = 0;
        std::shuffle(mask.data(), mask.data() + mask.size(), rng_inject_errs);
        
        // inject the errors where requested
        for(int i = 0; i < word.size(); i++)
            word[i] ^= mask[i];
    }
    else
    {
        assert(n_errors == 0 && "impossible to meet requests for the given configuration");
        word = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(word.size());
        return;
    }
}

/**
 * @brief convert rapidjson object to a printable C++-style string
 * 
 * @tparam T RapidJson type to stringify
 * @param o object to stringify
 * @return std::string C++-style string representation of object `o'
 */
template<typename T>
std::string rapidjson_stringify(const T& o)
{
    rapidjson::StringBuffer sb;
    rapidjson::Writer< rapidjson::StringBuffer > writer(sb);
    o.Accept(writer);
    return sb.GetString();
}

std::vector< std::vector< einsim::error_model_descriptor * > > einsim::error_model_descriptors_from_json(const std::string &json_cfg_file)
{

    // read the configuration file and determine how to construct this object
    std::ifstream ifs(json_cfg_file);
    rapidjson::IStreamWrapper isw(ifs);

    rapidjson::Document d;
    rapidjson::ParseResult result = d.ParseStream< rapidjson::kParseCommentsFlag >(isw);
    if(!result) 
    {
        std::cout << "[ERROR] JSON parse error while parsing configuration file " <<
            json_cfg_file << ": " << rapidjson::GetParseError_En(result.Code()) 
            << " (" << result.Offset() << ")" << std::endl;
        assert(0 && "malformed JSON configuration file");
    }

    // extract each of the error models
    if(!d.IsArray())
    {
        std::cout << "[ERROR] malformed error model configuration file " <<
            json_cfg_file << ": top level must be an array" << std::endl;
        assert(0 && "malformed JSON configuration file");
    }

    std::vector< std::vector< einsim::error_model_descriptor * > > error_models;
    for(const auto &emd_vec_cfg : d.GetArray())
    {
        std::vector< einsim::error_model_descriptor * > emd_vec;
    
        if(!emd_vec_cfg.IsArray())
        {
            std::cout << "[ERROR] malformed error model configuration file " <<
                json_cfg_file << ": second level must be an array" << std::endl;
            assert(0 && "malformed JSON configuration file");
        }

        // create EMD's per bit
        std::vector< std::vector< einsim::error_model_descriptor * > > emds_per_bit;
        for(const auto &bit_emd_cfg : emd_vec_cfg.GetArray())
        {
            const std::string &bit_error_model = bit_emd_cfg["error_model"].GetString();
            const auto &bit_model_params_array = bit_emd_cfg["model_params"].GetArray();
            std::vector< einsim::error_model_descriptor * > emds_this_bit;
            einsim::error_model em = einsim::str_to_enum_error_model(bit_error_model);
            for(const auto &param_set : bit_model_params_array)
            {
                // create an EMD for this param set
                std::vector< std::string > model_params;
                for(const auto &param : param_set.GetArray())
                    model_params.push_back(rapidjson_stringify(param));
                einsim::error_model_descriptor *emd = einsim::error_model_descriptor_from_params(em, model_params);
                emds_this_bit.push_back(emd);
            }
            emds_per_bit.push_back(emds_this_bit);
        }

        // now, form the cartesian product of all bit sets and and it to the final set
        einsim::construct_cartesian_product_of_per_bit_error_models(error_models, emds_per_bit);
    }
    return error_models;
}

void einsim::construct_cartesian_product_of_per_bit_error_models(
      std::vector< std::vector< einsim::error_model_descriptor * > > &error_models
    , const std::vector< std::vector< einsim::error_model_descriptor * > > &emds_per_bit)
{
    std::cout << "[DEBUG] forming cartesian product of [";
    for(const auto &elem : emds_per_bit)
        std::cout << elem.size() << " ";
    std::cout << "]" << std::endl;

    // counters for each bit
    std::vector< size_t > count;
    for(size_t bit = 0; bit < emds_per_bit.size(); bit++)
        count.push_back(0);

    // go through each counter value until we hit the end
    int n = 0;
    while(true)
    {
        // product of the counters
        std::vector< einsim::error_model_descriptor * > bitvector;
        for(size_t bit = 0; bit < emds_per_bit.size(); bit++)
            bitvector.push_back(emds_per_bit.at(bit).at(count.at(bit)));
        error_models.push_back(bitvector);
        n++;

        // increment the counters
        assert(count.size() > 0);
        count.at(0)++;
        for(size_t bit = 0; bit < emds_per_bit.size() - 1; bit++)
        {
            if(count.at(bit) == emds_per_bit.at(bit).size())
            {
                count.at(bit) = 0;
                count.at(bit + 1)++;
            }
        }
        if(count.at(count.size() - 1) == emds_per_bit.at(emds_per_bit.size() - 1).size())
            break;
    }
    std::cout << "[DEBUG] resulted in " << n << " entries" << std::endl;
}

/**
 * @brief convert string representation of a boolean to the C++ bool type
 * 
 * @param str string representation of boolean
 * @return true true
 * @return false false
 */
bool stob(const std::string &str) 
{
    bool b;
    std::string str_lower;
    std::transform(str.begin(), str.end(), str_lower.begin(), ::tolower);
    std::istringstream is(str);
    is >> std::boolalpha >> b;
    return b;
}

einsim::error_model_descriptor *einsim::error_model_descriptor_from_params(enum error_model em, const std::vector< std::string > &model_params)
{
    std::cout << "[DEBUG] instantiating error model " << einsim::enum_to_str_error_model(em) 
        << " with params [";
    for(const std::string &s : model_params)
        std::cout << s << " ";
    std::cout << "]" << std::endl;
    
    switch(em)
    {
        case EM_NORMAL:
        {
            assert(model_params.size() == 0 && "invalid number of parameters for EM_NORMAL");
            einsim::error_model_descriptor *emd = new einsim::error_model_descriptor_normal();
            return emd;
        }

        case EM_UNIFORM_RANDOM:
        {
            assert(model_params.size() == 1 && "invalid number of parameters for EM_UNIFORM_RANDOM");
            float error_rate = std::stod(model_params.at(0));
            einsim::error_model_descriptor *emd = new einsim::error_model_descriptor_uniform_random(error_rate);
            return emd;
        }

        case EM_DATA_RETENTION:
        {
            assert(model_params.size() == 1 && "invalid number of parameters for EM_DATA_RETENTION");
            float error_rate = std::stod(model_params.at(0));
            einsim::error_model_descriptor *emd = new einsim::error_model_descriptor_data_retention(error_rate);
            return emd;
        }

        case EM_DATA_RETENTION_NOISY:
        {
            assert(model_params.size() == 2 && "invalid number of parameters for EM_DATA_RETENTION_NOISY");
            float error_rate = std::stod(model_params.at(0));
            float noise_ratio = std::stod(model_params.at(1));
            einsim::error_model_descriptor *emd = new einsim::error_model_descriptor_data_retention_noisy(error_rate, noise_ratio);
            return emd;
        }

        case EM_STUCK_AT:
        {
            assert(model_params.size() == 1 && "invalid number of parameters for EM_STUCK_AT");
            bool stuck_at_val = stob(model_params.at(0));
            einsim::error_model_descriptor *emd = new einsim::error_model_descriptor_stuck_at(stuck_at_val);
            return emd;
        }

        default:
        case EM_UNKNOWN:
            std::cout << "[ERROR] cannot construct EM_UNKNOWN error model from any params" << std::endl;
            assert(0);
    }   
}

int einsim::error_model_descriptor::get_n_model_params(enum error_model em)
{
    switch(em)
    {
        case EM_NORMAL:               return error_model_descriptor_normal::get_n_model_params();
        case EM_UNIFORM_RANDOM:       return error_model_descriptor_uniform_random::get_n_model_params();
        case EM_DATA_RETENTION:       return error_model_descriptor_data_retention::get_n_model_params();
        case EM_DATA_RETENTION_NOISY: return error_model_descriptor_data_retention_noisy::get_n_model_params();
        case EM_STUCK_AT:             return error_model_descriptor_stuck_at::get_n_model_params();
        default:
        case EM_UNKNOWN:
            std::cout << "[ERROR] cannot ask for number of params for EM_UNKNOWN error model" << std::endl;
            return -1;

    }
}

std::string einsim::error_model_descriptor_vec_to_str(const std::vector< einsim::error_model_descriptor * > &emd_vec)
{
    std::string str = "";
    for(const einsim::error_model_descriptor *emd : emd_vec)
        str += emd->to_str() + ";";
    return str.substr(0, str.size() - 1);
}
