/**
 * @file ecc_code.cpp
 *
 * @brief Supporting routines for the ecc_code class
 *
 * @author Minesh Patel
 * Contact: minesh.patelh@gmail.com
 */
#include <iostream>
#include <string>
#include <algorithm>
#include <random>
#include <fstream>

/* libraries */
#include "libtp.h"
#include "Eigen/Eigen"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"

/* project files */
#include "word_generator.h"
#include "error_model.h"
#include "supporting_routines.h"
#include "ecc_code.h"
#include "codes/repetition_code.h"
#include "codes/bch_code.h"
#include "codes/hamming_code.h"

std::string einsim::enum_to_str_test_mode(enum einsim::test_mode tm)
{
    static std::map< enum einsim::test_mode, std::string > enum_to_str_test_mode_map = 
    {
          {einsim::TM_FAST, "FAST"}
        , {einsim::TM_SLOW, "SLOW"}
    };
    auto search = enum_to_str_test_mode_map.find(tm);
    if(search == enum_to_str_test_mode_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

enum einsim::test_mode einsim::str_to_enum_test_mode(const std::string &str)
{
    static std::map< std::string, enum einsim::test_mode > str_to_enum_test_mode_map = 
    {
          {"FAST", einsim::TM_FAST}
        , {"SLOW", einsim::TM_SLOW}
    };
    std::string str_uppercase(str);
    std::transform(str_uppercase.begin(), str_uppercase.end(), str_uppercase.begin(), ::toupper);
    auto search = str_to_enum_test_mode_map.find(str_uppercase);
    if(search == str_to_enum_test_mode_map.end())
        return einsim::TM_UNKNOWN;
    else
        return search->second;
}

std::string einsim::get_all_possible_test_modes(void)
{
    std::string ret;
    for(int m = einsim::TM_FAST; m < einsim::TM_UNKNOWN; m++)
        ret += einsim::enum_to_str_test_mode((einsim::test_mode)m) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}

std::string einsim::enum_to_str_ecc_scheme(enum einsim::ecc_scheme es)
{
    static std::map< enum einsim::ecc_scheme, std::string > enum_to_str_ecc_scheme_map = 
    {
          {einsim::ES_REPETITION_T1, "REP_T1"}
        , {einsim::ES_REPETITION_T2, "REP_T2"}
        , {einsim::ES_REPETITION_T3, "REP_T3"}
        , {einsim::ES_HAMMING_SEC, "HSC"}
        , {einsim::ES_BCH_T1, "BCH_T1"}
        , {einsim::ES_BCH_T2, "BCH_T2"}
        , {einsim::ES_BCH_T3, "BCH_T3"}
    };
    auto search = enum_to_str_ecc_scheme_map.find(es);
    if(search == enum_to_str_ecc_scheme_map.end())
        return "UNKNOWN";
    else
        return search->second;
}

enum einsim::ecc_scheme einsim::str_to_enum_ecc_scheme(const std::string &es)
{
    static std::map< std::string, enum ecc_scheme > str_to_enum_ecc_scheme_map = 
    {
          {"REP_T1", einsim::ES_REPETITION_T1}
        , {"REP_T2", einsim::ES_REPETITION_T2}
        , {"REP_T3", einsim::ES_REPETITION_T3}
        , {"HSC",    einsim::ES_HAMMING_SEC}
        , {"BCH_T1", einsim::ES_BCH_T1}
        , {"BCH_T2", einsim::ES_BCH_T2}
        , {"BCH_T3", einsim::ES_BCH_T3}
    };
    std::string es_uppercase(es);
    std::transform(es_uppercase.begin(), es_uppercase.end(), es_uppercase.begin(), ::toupper);
    auto search = str_to_enum_ecc_scheme_map.find(es_uppercase);
    if(search == str_to_enum_ecc_scheme_map.end())
    {
        return einsim::ES_UNKNOWN;
    }
    else
        return search->second;

}

std::string einsim::get_all_possible_ecc_schemes(void)
{
    std::string ret;
    for(int s = einsim::ES_REPETITION_T1; s < einsim::ES_UNKNOWN; s++)
        ret += einsim::enum_to_str_ecc_scheme((einsim::ecc_scheme)s) + ", ";
    ret.erase(ret.size() - 2, ret.size() - 1);
    return ret;
}

einsim::ecc_code *einsim::build_ecc_code(enum einsim::ecc_scheme scheme, int n_data_bits, int random_seed)
{
    einsim::ecc_code *ecc_scheme_p = NULL;
    switch(scheme)
    {
        case einsim::ES_REPETITION_T1:
        case einsim::ES_REPETITION_T2:
        case einsim::ES_REPETITION_T3:
        {
            int correction_capability = (scheme == einsim::ES_REPETITION_T1) ? 3 : ((scheme == einsim::ES_REPETITION_T2) ? 5 : 7);
            ecc_scheme_p = new einsim::repetition(random_seed, n_data_bits, correction_capability);
            break;
        }
        case einsim::ES_HAMMING_SEC:
        {
            ecc_scheme_p = new einsim::hamming(random_seed, n_data_bits);
            break;
        }
        case einsim::ES_BCH_T1:
        case einsim::ES_BCH_T2:
        case einsim::ES_BCH_T3:
        {
            int correction_capability = (scheme == einsim::ES_BCH_T1) ? 1 : ((scheme == einsim::ES_BCH_T2) ? 2 : 3);
            einsim::bch *bch_code = new einsim::bch(random_seed, n_data_bits, correction_capability);
            if(bch_code->ready())
                ecc_scheme_p = bch_code;
            else
            {
                printf_both("[ERROR] no such BCH code exists for parameters k: %d, t: %d!\n", n_data_bits, correction_capability);
                assert(0);
            }
            break;
        }
        default:
            printf("[ERROR] unknown/invalid ECC scheme requested: %s\n", enum_to_str_ecc_scheme(scheme).c_str());
            return NULL;            
    }

    return ecc_scheme_p;
}

einsim::ecc_code *einsim::build_ecc_code(const std::string &cfg_file_name)
{

    std::ifstream ifs(cfg_file_name);
    rapidjson::IStreamWrapper isw(ifs);

    rapidjson::Document d;
    rapidjson::ParseResult result = d.ParseStream< rapidjson::kParseCommentsFlag >(isw);
    if(!result) 
    {
        std::cout << "[ERROR] JSON parse error while parsing configuration file " <<
            cfg_file_name << ": " << rapidjson::GetParseError_En(result.Code()) 
            << " (" << result.Offset() << ")" << std::endl;
        return NULL;
    }

    // get the ECC type
    const std::string &ecc_scheme_str = d["s"].GetString();
    enum einsim::ecc_scheme ecc_scheme = str_to_enum_ecc_scheme(ecc_scheme_str);
    switch(ecc_scheme)
    {
        case einsim::ES_REPETITION_T1:
        case einsim::ES_REPETITION_T2:
        case einsim::ES_REPETITION_T3:
            assert(0 && "building repetition code from cfg file unsupported");

        case einsim::ES_HAMMING_SEC:
            return new einsim::hamming(d, cfg_file_name);

        case einsim::ES_BCH_T1:
        case einsim::ES_BCH_T2:
        case einsim::ES_BCH_T3:
            assert(0 && "building BCH code from cfg file unsupported");

        default:
            printf("[ERROR] unknown/invalid ECC scheme requested: %s\n", ecc_scheme_str.c_str());
            return NULL;            
    }

    // // pretty-print the JSON
    // rapidjson::StringBuffer buffer;
    // rapidjson::Writer< rapidjson::StringBuffer > writer(buffer);
    // d.Accept(writer);
    // std::cout << buffer.GetString() << std::endl;

    return NULL;
}


void einsim::test_ecc(void (*test_func)(thread_pool &tp, enum test_mode mode), enum test_mode mode, int n_threads)
{
    thread_pool tp(n_threads);
    tp.start();

    test_func(tp, mode);

    // wait for jobs to complete
    while(tp.get_n_jobs_outstanding())
    {   
        int n_total_jobs = tp.get_n_jobs_outstanding() + tp.get_n_jobs_completed();
        fprintf(g_output_file, "Testing: [%d/%d] jobs remaining\n", tp.get_n_jobs_outstanding(), n_total_jobs);
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
    tp.wait();
    
    printf_both("\nTest complete\n");
}

void einsim::test_thread(int tid, einsim::ecc_code *ec)
{
    // Generate a random dataword
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_sent = 
        Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */>::Zero(ec->get_n_data_bits());
    Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 /* cols */> dummy;
    einsim::generate_word(data_word_sent, DP_ALL_ONES, dummy, einsim::CD_ALL_TRUE);

    // Encode the dataword into codeword
    Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word_sent = ec->encode(data_word_sent);

    // Corrupt codeword with N errors
    for(int nerrs_transmission = 0; nerrs_transmission <= ec->get_n_code_bits(); nerrs_transmission++)
    {
        bool fully_correctable = ec->correction_capability() >= nerrs_transmission;
        
        // inject precisely N errors into the codeword (i.e., errors during transmission)
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word_recv = code_word_sent;
        inject_n(code_word_recv, einsim::EM_UNIFORM_RANDOM, einsim::CD_ALL_TRUE, einsim::DP_ALL_ONES, nerrs_transmission);

        // sanity-check that this is the right amount of errors
        int nerrs_induced = hamming_distance(code_word_sent, code_word_recv);
        if(nerrs_induced != nerrs_transmission)
        {
            std::stringstream ss;
            ss << "code_sent: " << code_word_sent.transpose() << std::endl;
            ss << "code_rcvd: " << code_word_recv.transpose() << std::endl;
            ss << "data_sent: " << data_word_sent.transpose() << std::endl;
            printf_both("[ERROR] #errs induced (%d) != #errs transmitted(%d) (%s)\n%s", nerrs_induced, nerrs_transmission, ec->name().c_str(), ss.str().c_str());
            exit(-1);
        }

        // Decode codeword into dataword
        Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word_recv = ec->decode(code_word_recv);

        // std::cout << "Dataword orig: " << data_word_sent.transpose() << std::endl;
        // std::cout << "Codeword sent: " << code_word_sent.transpose() << std::endl;
        // std::cout << "Codeword recv: " << code_word_recv.transpose() << std::endl;
        // std::cout << "Dataword recv: " << data_word_recv.transpose() << std::endl;

        int nerrs_after = hamming_distance(data_word_sent, data_word_recv);
        if(fully_correctable && nerrs_after)
        {
            std::stringstream ss;
            ss << "    > code_sent: " << code_word_sent.transpose() << std::endl;
            ss << "    > code_rcvd: " << code_word_recv.transpose() << std::endl;
            ss << "    > data_sent: " << data_word_sent.transpose() << std::endl;
            ss << "    > data_rcvd: " << data_word_recv.transpose() << std::endl;
            printf_both("[ERROR]: observed %d errors when %d induced and %d correctable (%s)\n%s", nerrs_after, nerrs_transmission, ec->correction_capability(), ec->name().c_str(), ss.str().c_str());
            exit(-1);
        }
    
        // if(nerrs_after)
        //     printf("nerrs: %d -> %d\n", nerrs_transmission, nerrs_after);
    }
}