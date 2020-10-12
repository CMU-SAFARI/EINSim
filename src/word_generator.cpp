/**
 * @file word_generator.cpp
 *
 * @brief Routines for generating words according to specific data patterns
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <string>
#include <algorithm>
#include <random>

/* libraries */
#include "Eigen/Eigen"

/* project includes */
#include "supporting_routines.h"
#include "word_generator.h"

static std::default_random_engine true_anti_cell_generator_rng;
static std::default_random_engine word_generator_rng;
static std::bernoulli_distribution bernoulli_0_5(0.5);
static bool word_generator_initialized = false;

enum einsim::data_pattern einsim::str_to_enum_data_pattern(const std::string &str)
{
    static std::map< std::string, enum data_pattern > str_to_enum_data_pattern_map = 
    {
          {"RANDOM", DP_RANDOM}
        , {"ALL_ONES", DP_ALL_ONES}
        , {"CHARGED", DP_CHARGED}
        , {"CUSTOM", DP_CUSTOM}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_data_pattern_map.find(str_uppercase);
    if(search == str_to_enum_data_pattern_map.end())
    {
        if(    !strncmp(str_uppercase.c_str(), "0B", 2)
            || !strncmp(str_uppercase.c_str(), "0O", 2)
            || !strncmp(str_uppercase.c_str(), "0X", 2))
            return DP_CUSTOM;
        else
            return DP_UNKNOWN;
    }
    else
        return search->second;
}

std::string einsim::enum_to_str_data_pattern(enum einsim::data_pattern dp)
{
    static std::map< enum data_pattern, std::string > enum_to_str_data_pattern_map = 
    {
          {DP_RANDOM, "RANDOM"}
        , {DP_ALL_ONES, "ALL_ONES"}
        , {DP_CHARGED, "CHARGED"}
        , {DP_CUSTOM, "CUSTOM"}
    };
    auto search = enum_to_str_data_pattern_map.find(dp);
    if(search == enum_to_str_data_pattern_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::custom_dp_to_vector(const std::string &str)
{
    // std::cout << "[DEBUG] custom_dp_to_vector input: " << str << std::endl;

    int base_log2 = 0;
    int base = 0;
    switch(str.at(1))
    {
        case 'b': case 'B': base = 2; base_log2 = 1; break;
        case 'o': case 'O': base = 8; base_log2 = 3; break;
        case 'x': case 'X': base = 16; base_log2 = 4; break;
        default:
        {                
            printf("[ERROR] Invalid custom data pattern: %s\n", str.c_str());
            assert(0 && "illegal data pattern - must start with 0b, 0o, 0x");
        }
    }

    // fill in the binary vector as we go, LSB onwards
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > dp;
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > symbol = Eigen::Matrix< ET, Eigen::Dynamic, 1 >::Zero(base_log2);
    for(int i = str.size() - 1; i >= 2; i--)
    {
        const char raw = str.at(i);
        std::string val_str = "";
        val_str += raw;
        char val = std::stoi(val_str, 0, base); // TODO: optimize using LUT or somesuch, but it only happens once
        for(int j = 0; j < base_log2; j++)
            symbol[base_log2 - 1 - j] = (val >> j) & 1;
        
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > new_dp;
        new_dp.resize(dp.rows() + symbol.rows(), 1);
        new_dp << symbol, dp;
        dp = new_dp;
        // std::cout << "DP: " << dp << " STR: " << str  << " VAL: " << val << " VAL_STR: " << val_str << std::endl;
    }

    // std::cout << "[DEBUG] custom_dp_to_vector output: " << dp << std::endl;
    return dp;
}

std::string einsim::custom_dp_to_str(Eigen::Matrix< ET, Eigen::Dynamic, 1 > dp)
{
    static char hex_lut[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};
    std::string ret;
    int vec_len = dp.rows();
    uint8_t x_val = 0;
    for(int bit_idx = 0; bit_idx < vec_len; bit_idx++)
    {
        x_val += dp[vec_len - bit_idx - 1] * (1 << (bit_idx % 4));
        if(!((bit_idx + 1) % 4) || (bit_idx == vec_len - 1))
        {
            ret = std::string(1, hex_lut[x_val]) + ret;
            x_val = 0;
        }
    }
    return ret;
}

std::string einsim::get_all_possible_data_patterns(void)
{
    std::string ret;
    for(int s = einsim::DP_RANDOM; s < einsim::DP_UNKNOWN; s++)
        ret += einsim::enum_to_str_data_pattern((einsim::data_pattern)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}

enum einsim::true_anti_cell_distribution einsim::str_to_enum_true_anti_cell_distribution(const std::string &str)
{
    static std::map< std::string, enum true_anti_cell_distribution > str_to_enum_true_anti_cell_distribution_map = 
    {
          {"ALL_TRUE_OR_ALL_ANTI", CD_ALL_TRUE_OR_ALL_ANTI}
        , {"ALL_TRUE", CD_ALL_TRUE}
        , {"ALL_ANTI", CD_ALL_ANTI}
        , {"COLSTRIPE_T", CD_COLSTRIPE_T}
        , {"COLSTRIPE_A", CD_COLSTRIPE_A}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_true_anti_cell_distribution_map.find(str_uppercase);
    if(search == str_to_enum_true_anti_cell_distribution_map.end())
        return CD_UNKNOWN;
    else
        return search->second;
}

std::string einsim::enum_to_str_true_anti_cell_distribution(enum einsim::true_anti_cell_distribution cd)
{
    static std::map< enum true_anti_cell_distribution, std::string > enum_to_str_true_anti_cell_distribution_map = 
    {
          {CD_ALL_TRUE_OR_ALL_ANTI, "ALL_TRUE_OR_ALL_ANTI"}
        , {CD_ALL_TRUE, "ALL_TRUE"}
        , {CD_ALL_ANTI, "ALL_ANTI"}
        , {CD_COLSTRIPE_T, "COLSTRIPE_T"}
        , {CD_COLSTRIPE_A, "COLSTRIPE_A"}
    };
    auto search = enum_to_str_true_anti_cell_distribution_map.find(cd);
    if(search == enum_to_str_true_anti_cell_distribution_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

std::string einsim::get_all_possible_true_anti_cell_distributions(void)
{
    std::string ret;
    for(int s = einsim::CD_ALL_TRUE_OR_ALL_ANTI; s < einsim::CD_UNKNOWN; s++)
        ret += einsim::enum_to_str_true_anti_cell_distribution((einsim::true_anti_cell_distribution)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}

enum einsim::true_anti_cell_state einsim::str_to_enum_true_anti_cell_state(const std::string &str)
{
    static std::map< std::string, enum true_anti_cell_state > str_to_enum_true_anti_cell_state_map = 
    {
          {"ALL_TRUE", TACS_ALL_TRUE}
        , {"ALL_ANTI", TACS_ALL_ANTI}
        , {"ALT_T",    TACS_ALT_T}
        , {"ALT_A",    TACS_ALT_A}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_true_anti_cell_state_map.find(str_uppercase);
    if(search == str_to_enum_true_anti_cell_state_map.end())
        return TACS_UNKNOWN;
    else
        return search->second;
}

std::string einsim::enum_to_str_true_anti_cell_state(enum einsim::true_anti_cell_state cd)
{
    static std::map< enum true_anti_cell_state, std::string > enum_to_str_true_anti_cell_state_map = 
    {
          {TACS_ALL_TRUE, "ALL_TRUE"}
        , {TACS_ALL_ANTI, "ALL_ANTI"}
        , {TACS_ALT_T,    "ALT_T"}
        , {TACS_ALT_A,    "ALT_A"}
    };
    auto search = enum_to_str_true_anti_cell_state_map.find(cd);
    if(search == enum_to_str_true_anti_cell_state_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

std::string einsim::get_all_possible_true_anti_cell_states(void)
{
    std::string ret;
    for(int s = einsim::TACS_ALL_TRUE; s < einsim::TACS_UNKNOWN; s++)
        ret += einsim::enum_to_str_true_anti_cell_state((einsim::true_anti_cell_state)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}


enum einsim::word_to_burst_mapping einsim::str_to_enum_word_to_burst_mapping(const std::string &str)
{
    static std::map< std::string, enum word_to_burst_mapping > str_to_enum_word_to_burst_mapping_map = 
    {
          {"BLOCKS", W2BM_BLOCKS}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_word_to_burst_mapping_map.find(str_uppercase);
    if(search == str_to_enum_word_to_burst_mapping_map.end())
        return W2BM_UNKNOWN;
    else
        return search->second;
}

std::string einsim::enum_to_str_word_to_burst_mapping(enum einsim::word_to_burst_mapping cd)
{
    static std::map< enum word_to_burst_mapping, std::string > enum_to_str_word_to_burst_mapping_map = 
    {
          {W2BM_BLOCKS, "BLOCKS"}
    };
    auto search = enum_to_str_word_to_burst_mapping_map.find(cd);
    if(search == enum_to_str_word_to_burst_mapping_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

std::string einsim::get_all_possible_word_to_burst_mappings(void)
{
    std::string ret;
    for(int s = einsim::W2BM_BLOCKS; s < einsim::W2BM_UNKNOWN; s++)
        ret += einsim::enum_to_str_word_to_burst_mapping((einsim::word_to_burst_mapping)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}

einsim::true_anti_cell_state einsim::generate_word(
      Eigen::Matrix< ET, Eigen::Dynamic, 1 > &word
    , const enum einsim::data_pattern dp
    , const Eigen::Matrix< ET, Eigen::Dynamic, 1 > &custom_dp
    , const enum einsim::true_anti_cell_distribution cd)
{
    // handle the one-time initialization of the global random number generator
    if(!word_generator_initialized)
    {
        true_anti_cell_generator_rng.seed(std::random_device{}());
        word_generator_rng.seed(std::random_device{}());
        word_generator_initialized = true;
    }

    // determine the true-/anti-cell state of this word
    einsim::true_anti_cell_state tacs;
    switch(cd)
    {
        case CD_ALL_TRUE_OR_ALL_ANTI:
            tacs = bernoulli_0_5(true_anti_cell_generator_rng) ? TACS_ALL_TRUE : TACS_ALL_ANTI;
            break;

        case CD_ALL_TRUE:
            tacs = TACS_ALL_TRUE;
            break;

        case CD_ALL_ANTI:
            tacs = TACS_ALL_ANTI;
            break;

        case CD_COLSTRIPE_T:
            tacs = TACS_ALT_T;
            break;

        case CD_COLSTRIPE_A:
            tacs = TACS_ALT_A;
            break;

        default:
            printf_both("[ERROR] Unsupported true-/anti-cell distribution: %d\n", (int)cd);
            assert(0 && "unsupported true-/anti-cell distribution");
    }

    // write in the appropriate data pattern to the word
    // may or may not be sensitive to the true-/anti-cell layout, so be careful when adding new entries
    switch(dp)
    {
        case einsim::DP_RANDOM:
        {
            for(int i = 0; i < word.size(); i++)
                word[i] = bernoulli_0_5(word_generator_rng) ? 0 : 1;
            break;
        }
        case einsim::DP_ALL_ONES:
        {
            word = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Ones(word.size());
            break;
        }
        case einsim::DP_CHARGED:
        {
            switch(tacs)
            {
                case einsim::TACS_ALL_TRUE:
                    word = Eigen::Matrix< ET, Eigen::Dynamic, 1>::Ones(word.rows());
                    break;

                case einsim::TACS_ALL_ANTI:
                    word = Eigen::Matrix< ET, Eigen::Dynamic, 1>::Zero(word.rows());
                    break;
                    
                case einsim::TACS_ALT_T:
                    word = Eigen::Matrix< ET, Eigen::Dynamic, 1>::Zero(word.rows());
                    for(int i = 0; i < word.rows(); i += 2)
                        word(i) = 1;
                    break;
                    
                case einsim::TACS_ALT_A:
                    word = Eigen::Matrix< ET, Eigen::Dynamic, 1>::Ones(word.rows());
                    for(int i = 0; i < word.rows(); i += 2)
                        word(i) = 0;
                    break;
                    
                default:
                    printf_both("[ERROR] Unsupported true-/anti-cell state: %d\n", (int)tacs);
                    assert(0 && "unsupported true-/anti-cell state");
            }
            break;
        }
        case einsim::DP_CUSTOM:
        {
            assert(word.rows() == custom_dp.rows() && word.cols() == custom_dp.cols() && "internal error- custom_dp should match word dimensions");
            word = custom_dp;
            break;
        }
        default:
        {
            printf_both("[ERROR] unknown or unsupported data pattern: %d", dp);
            exit(-1);
        }
    }

    return tacs;
}
