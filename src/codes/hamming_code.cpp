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
#include <regex>

/* libraries */
#include "Eigen/Eigen"
#include "libtp.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"

/* project includes */
#include "supporting_routines.h"
#include "ecc_code.h"
#include "hamming_code.h"

/**
 * @brief compute the number of parity-check bits for a Hamming code
 * 
 * @param n_data_bits number of data bits in the code
 * @return int number of parity-check bits
 */
int compute_n_parity_bits(int n_data_bits)
{
    int n_parity_bits = 0;
    while((1 << n_parity_bits) < (n_parity_bits + n_data_bits + 1))
        n_parity_bits++;
    assert(n_parity_bits == (int) std::ceil(std::log2(n_data_bits + n_parity_bits + 1)));
    return n_parity_bits;
}

/**
 * @brief tests whether a positive integer is a power of 2
 * 
 * @param n positive integer to test
 * @return true n is a power of 2
 * @return false n is not a power of 2
 */
bool is_po2(int n)
{
    assert(n > 0 && "does not work for zero or negative numbers");
    return (n & (n - 1)) == 0;
}

/**
 * @brief computes the {G, H, R} matrices for a Hamming code given the input parameters
 * 
 * @param nd number of data bits
 * @param np number of parity-check bits
 * @param permutation random seed
 * @param G generator matrix
 * @param H parity-check matrix
 * @param R degenerator matrix
 * @param data_bit_indices list of bit positions corresponding to data bits
 * @param parity_bit_indices list of bit positions corresponding to parity-check bits
 * @param use_standard_form whether to keep the G/H/R matrices in standard form
 */
void compute_hamming_matrices
(
      int nd, int np, int permutation
    , Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> &G
    , Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> &H
    , Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> &R
    , std::vector< int > data_bit_indices
    , std::vector< int > parity_bit_indices
    , bool use_standard_form)
{
    // generate a set of each possibles syndrome value up to total number of bits
    std::vector< int > syn_values_non_po2;
    std::vector< int > syn_values_po2;
    for(int i = 1; i < (1 << np); i++)
        if(i & (i - 1))
            syn_values_non_po2.push_back(i);
        else
            syn_values_po2.push_back(i);
    assert(syn_values_po2.size() == (unsigned)np && "mismatch in computed number of po2 syndromes");

    // select a random subset of (nd) non-po2 items
    std::seed_seq seed{permutation};
    std::mt19937 rng(seed);
    std::vector< int > indices(syn_values_non_po2.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.data(), indices.data() + indices.size(), rng);

    // calculate the final syndrome values - form [I_{n-k} | P^t]
    std::vector< int > syn_values(syn_values_po2);
    for(int i = 0; i < nd; i++)
        syn_values.push_back(syn_values_non_po2[indices[i]]);

    // randomize the final syn values order
    indices.resize(syn_values.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.data(), indices.data() + indices.size(), rng);

    // construct the parity-check matrix based on the final syn value order
    H = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(np, nd + np);
    assert(H.rows() == np);
    assert(H.cols() == np + nd);         
    for(int col = 0; col < nd + np; col++) // inject syndromes
        for(int row = 0; row < np; row++)
            H(row, col) = (syn_values[indices[col]] >> row) & 1; 

    // identify data bit indices and parity bit indices in the codeword
    data_bit_indices.clear();
    parity_bit_indices.clear();
    for(int i = 0; i < nd + np; i++)
    {
        if(is_po2(syn_values.at(indices.at(i))))
            parity_bit_indices.push_back(i);
        else
            data_bit_indices.push_back(i);
    }
    
    // row-reduce the parity-check matrix to obtain the coordinate vectors in the columns
    Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> parity_check_rref = row_reduce_to_rref(H, 0);

    // determine the permutation matrix that rearranges the the columns into standard form [P^t | I_{n-k}]
    Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > permutation_matrix = 
        Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(nd + np, nd + np);
    for(int syn_idx = 0; syn_idx < np; syn_idx++)
    {
        int syn = (1 << syn_idx);
        Eigen::Matrix< ET, Eigen::Dynamic /* rows */, 1 > syn_as_vec(np);
        for(int i = 0; i < np; i++)
            syn_as_vec(i) = (syn >> i) & 1;
        for(int c = 0; c < nd + np; c++)
            if(parity_check_rref.col(c).isApprox(syn_as_vec))
            {
                swap_matrix_rows(permutation_matrix, c, nd + syn_idx);
                break;
            }
    }
    Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> parity_check_standard_form = parity_check_rref * permutation_matrix;

    // construct the generator matrix using the form [I_k | P]^t
    Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> generator_standard_form 
        = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(nd + np, nd);
    assert(generator_standard_form.rows() == nd + np);
    assert(generator_standard_form.cols() == nd);
    generator_standard_form.block(nd, 0, np, nd) = parity_check_standard_form.block(0, 0, np, nd);
    generator_standard_form.block(0, 0, nd, nd) = Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic >::Identity(nd, nd);
    
    // construct the degenerator matrix (rows represent order of bits to choose)
    Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> degenerator_standard_form =
        Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(nd, nd + np);
    degenerator_standard_form.block(0, 0, nd, nd) = Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic >::Identity(nd, nd);

    // apply the permutation order to G/R
    G = permutation_matrix * generator_standard_form;
    R = degenerator_standard_form * permutation_matrix.transpose();

    // force standard form
    if(use_standard_form)
    {
        H = parity_check_standard_form;
        G = generator_standard_form;
        R = degenerator_standard_form;
    }

    // print the G/H before permutation
    if(g_verbosity >= 3)
    {
        std::cout << "[INFO] Created Hamming parity matrices:" << std::endl;
        std::cout << "[INFO]     H:" << std::endl << H << std::endl;
        std::cout << "[INFO]     RREF(H):" << std::endl << parity_check_rref << std::endl;
        std::cout << "[INFO]     P:" << std::endl << permutation_matrix << std::endl;
        std::cout << "[INFO]     P * RREF(H):" << std::endl << parity_check_rref * permutation_matrix << std::endl;
        std::cout << "[INFO]     G:" << std::endl << generator_standard_form << std::endl;
        std::cout << "[INFO]     R:" << std::endl << degenerator_standard_form << std::endl;
        std::cout << "[INFO]     P * G:" << std::endl << G << std::endl;
        std::cout << "[INFO]     R * P:" << std::endl << R << std::endl;
    }
}

/**
 * @brief read a matrix out of a rapidjson object
 * 
 * @param json_obj rapidjson object representing a matrix
 * @param mat matrix to read out
 * @return int return code 0 on success, !0 on error
 */
int rapidjson_read_matrix(const rapidjson::Value &json_obj, Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */> &mat)
{
    // read the rows out of the file
    std::vector< Eigen::Matrix< ET, 1 /* rows */, Eigen::Dynamic /* cols */> > mat_rows;
    std::set< int > row_lengths;
    for(const auto &json_row : json_obj.GetArray())
    {
        Eigen::Matrix< ET, 1 /* rows */, Eigen::Dynamic /* cols */> row;
        for(const auto &json_element : json_row.GetArray())
        {
            row.conservativeResize(Eigen::NoChange, row.cols() + 1);
            row(row.cols() - 1) = json_element.GetInt();
        }
        mat_rows.push_back(row);
        row_lengths.insert(row.cols());
    }

    // ensure that all rows are of the same number of columns
    if(row_lengths.size() != 1)
    {
        std::cout << "[ERROR] matrix rows must be same the length. have:" << std::endl;
        for(const auto &e : row_lengths)
            std::cout << "    " << e << std::endl;
        return -1;
    }

    // build the matrix to return
    int ncols = *row_lengths.begin();
    int nrows = mat_rows.size();
    mat = Eigen::Matrix< ET, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(nrows, ncols);
    for(int r = 0; r < nrows; r++)
        mat.row(r) = mat_rows.at(r);
    return 0;
}

/**
 * @brief sanity-checks the {G, H} matrices against the basic properties of Hamming codes
 * 
 * @param G generator matrix
 * @param H parity-check matrix
 * @return int return code 0 on success, !0 on error
 */
int check_hamming_matrices(const Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &G, const Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > &H)
{
    // ensure G * H = 0
    // (n, k) * (p, n)
    Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > HG = MOD2_EIGEN_VEC(H * G);
    if(!HG.isZero((ET)1e-3)) // tolerance
    {
        std::cout << "[ERROR] invalid Hamming matrices: H*G != 0" << std::endl;
        std::cout << "H:" << std::endl << H << std::endl;
        std::cout << "G:" << std::endl << G << std::endl;
        std::cout << "HG:" << std::endl << HG << std::endl;
        return -1;
    }

    // ensure columns of H are nonzero
    for(int i = 0; i < H.cols(); i++)
        if(H.col(i).isZero((ET)1e-3))
        {
            std::cout << "[ERROR] invalid Hamming matrices: columns of H are NOT nonzero" << std::endl;
            std::cout << "H:" << std::endl << H << std::endl;
            std::cout << "G:" << std::endl << G << std::endl;
            std::cout << "HG:" << std::endl << HG << std::endl;       
            return -1;
        }

    // ensure columns of H are unique
    for(int i = 0; i < H.cols(); i++)
        for(int j = i + 1; j < H.cols(); j++)
            if(H.col(i).isApprox(H.col(j)))
            {
                std::cout << "[ERROR] invalid Hamming matrices: columns of H are NOT unique" << std::endl;
                std::cout << "H:" << std::endl << H << std::endl;
                std::cout << "G:" << std::endl << G << std::endl;
                std::cout << "HG:" << std::endl << HG << std::endl;       
                return -1;
            }

    return 0;
}

einsim::hamming::hamming(int permutation, int n_data_bits)
    : initialized(false), permutation(permutation), nd(n_data_bits) 
{
    bool use_standard_form = true;
    if(n_data_bits <= 0)
    {
        printf_both("ERROR: Invalid number of data bits: %d\n", n_data_bits);
        exit(-1);
    }

    // compute n_parity_bits
    np = compute_n_parity_bits(n_data_bits);
    std::cout << "[INFO] Generating hamming code of permutation " << permutation << " with: " << nd << " data bits, " << np << " parity bits" << std::endl;

    // compute the hamming matrices
    compute_hamming_matrices(nd, np, permutation, generator, parity_check, degenerator, data_bit_indices, parity_bit_indices, use_standard_form);

    // ensure that the code is created correctly
    if(check_hamming_matrices(generator, parity_check))
    {
        printf_both("[ERROR] invalid Hamming matrices\n");
        exit(-1);
    }

    // set the UID based on the matrices
    uid = compute_uid();

    initialized = true;
}

// initialize from configuration file
einsim::hamming::hamming(const rapidjson::Document &d, const std::string &cfg_file_name)
    : initialized(false) 
{
    assert(str_to_enum_ecc_scheme(d["s"].GetString()) == ES_HAMMING_SEC);
    permutation = d["p"].GetInt();
    nd = d["k"].GetInt();
    uid = d["uid"].GetUint64();

    if(nd <= 0)
    {
        printf_both("[ERROR] Invalid number of data bits: %d\n", nd);
        exit(-1);
    }

    // compute n_parity_bits
    np = compute_n_parity_bits(nd);
    std::cout << "[INFO] Reading hamming code of permutation " << permutation << " with: " << nd << " data bits, " << np << " parity bits from configuration file: " << cfg_file_name << std::endl;

    // read out the hamming matrices from the file
    if(d.HasMember("GT"))
    {
        Eigen::Matrix< ET, Eigen::Dynamic, Eigen::Dynamic > generator_T;
        if(rapidjson_read_matrix(d["GT"], generator_T))
        {
            std::cout << "[ERROR] unable to read GT matrix out of configuration file: " << cfg_file_name << std::endl;
            exit(-1);
        }
        generator = generator_T.transpose();
    }
    else if(rapidjson_read_matrix(d["G"], generator))
    {
        std::cout << "[ERROR] unable to read G matrix out of configuration file: " << cfg_file_name << std::endl;
        exit(-1);
    }
    
    if(rapidjson_read_matrix(d["H"], parity_check) || rapidjson_read_matrix(d["R"], degenerator))
    {
        std::cout << "[ERROR] unable to H/R read matrices out of configuration file: " << cfg_file_name << std::endl;
        exit(-1);
    }

    // ensure that the code is created correctly
    if(check_hamming_matrices(generator, parity_check))
    {
        printf_both("[ERROR] invalid Hamming matrices\n");
        exit(-1);
    }

    // check UID
    uint64_t computed_uid = compute_uid();
    if(computed_uid != uid)
    {
        printf_both("[ERROR] UID mismatch in %s code. cfg_file: %" PRIu64 " computed: %" PRIu64 "\n", static_name().c_str(), uid, computed_uid);
        assert(0 && "UID mismatch!");
    }

    initialized = true;
}

int einsim::hamming::to_json(std::string &json) const
{
    rapidjson::Document d;
    d.SetObject();
    rapidjson::Document::AllocatorType &allocator = d.GetAllocator();

    // populate the rapidjson document with necessary fields
    d.AddMember("s", "HSC", allocator);
    d.AddMember("k", get_n_data_bits(), allocator);
    d.AddMember("p", get_permutation(), allocator);
    d.AddMember("uid", (uint64_t)get_uid(), allocator);

    // GHR matrices
    rapidjson::Value G_mat(rapidjson::kArrayType);
    rapidjson::Value H_mat(rapidjson::kArrayType);
    rapidjson::Value R_mat(rapidjson::kArrayType);
    for(int r = 0; r < generator.rows(); r++)
    {
        rapidjson::Value G_mat_row(rapidjson::kArrayType);
        for(int c = 0; c < generator.cols(); c++)
            G_mat_row.PushBack(generator(r, c), allocator);
        G_mat.PushBack(G_mat_row, allocator);
    }
    for(int r = 0; r < parity_check.rows(); r++)
    {
        rapidjson::Value H_mat_row(rapidjson::kArrayType);
        for(int c = 0; c < parity_check.cols(); c++)
            H_mat_row.PushBack(parity_check(r, c), allocator);
        H_mat.PushBack(H_mat_row, allocator);
    }
    for(int r = 0; r < degenerator.rows(); r++)
    {
        rapidjson::Value R_mat_row(rapidjson::kArrayType);
        for(int c = 0; c < degenerator.cols(); c++)
            R_mat_row.PushBack(degenerator(r, c), allocator);
        R_mat.PushBack(R_mat_row, allocator);
    }
    d.AddMember("G", G_mat, allocator);
    d.AddMember("H", H_mat, allocator);
    d.AddMember("R", R_mat, allocator);

    // write out the string
    rapidjson::StringBuffer strbuf;
    rapidjson::PrettyWriter< rapidjson::StringBuffer > writer(strbuf);
    writer.SetFormatOptions(rapidjson::kFormatSingleLineArray);
    d.Accept(writer);
    
    json.clear();
    json = strbuf.GetString();
    
    // this regex_replace DOES NOT WORK on OSX for some reason, so instead we manually iterate
    // json = std::regex_replace(json, std::regex("], \\["), "]\n        , ["); // format arrays nicely
    size_t index = 0;
    while(true) 
    {
        index = json.find("], [", index);
        if(index == std::string::npos) 
            break;

        json.replace(index, 4, "]\n        , [");
        index += 13;
    }
    // std::cout << "JSON: " << json << std::endl;

    return 0;
}

uint64_t einsim::hamming::compute_uid(void) const
{
    return hash_matrix({generator, parity_check, degenerator});
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::hamming::encode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > data_word) const
{
    assert(data_word.size() == nd);
    return MOD2_EIGEN_VEC(generator * data_word);
}

Eigen::Matrix< ET, Eigen::Dynamic, 1 > einsim::hamming::decode(Eigen::Matrix< ET, Eigen::Dynamic, 1 > code_word) const
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
