#include "sieve.h"
#include "num_theory.h"
#include "util.h"

#include<gmp.h>
#include<gmpxx.h>

#include<stdint.h>
#include<unistd.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cassert>
#include<thread>
#include<string>
#include<cstdlib>

typedef unsigned long long ull;

typedef struct params_s {
    uint64_t B;
    mpz_class n;
    mpz_class chunk_size;
    mpz_class range;
    uint64_t threads;
} params_t;

bool parse_params(int argc, char *argv[], params_t &params) {
    int c;
    while ((c = getopt(argc, argv, "B:n:r:t:"))){
        switch(c){
            case 'B':
                params.B = std::stoi(optarg);
                break;
            case 'n':
                params.n = optarg;
                break;
            case 'r':
                params.range = optarg;
                break;
            case 't':
                params.threads = std::stoi(optarg);
                break;
        }
    }
    return true;
}

int main(int argc, char *argv[]) {
    //param_t params = {0,0,0,0,0};
    mpz_class n, sqrt_n;
    // RSA 129
    // n = "114381625757888867669235779976146612010218296721242362562561842935706935245733897830597123563958705058989075147599290026879543541";

    // RSA 100
    n = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
    sqrt_n = sqrt(n);

    // there are approx B/(2*log B) primes in the factor base
    uint64_t B = 40000000; // 10^7
    uint64_t chunk_size = 10*B;
    mpz_class range = 1000000000000; // 10^12

    std::cerr << "[INFO] B=" << B << std::endl;

    std::cerr << "[INFO] " << "n is " << mpz_sizeinbase(n.get_mpz_t(), 2)
              << " bits" << std::endl;

    std::vector<int64_t> primes = prime_sieve(B);
    std::vector<int64_t> factor_base;
    std::vector<int64_t> sqrts;

    for(const auto p: primes){
        int64_t n_modp = mpz_to_int64(n%p);
        if(p>=5 && is_residue(n_modp, p)){
            int64_t root = mod_sqrt(n_modp, p);
            factor_base.push_back(p);
            sqrts.push_back(root);
        }
    }
    std::cerr << "[INFO] " << factor_base.size() << " primes found" << std::endl;

    // sieving
    // std::ofstream outfile;
    // outfile.open("out.txt");
    // if (!outfile.is_open()) {
    //     std::cerr << "error opening out.txt" << std::endl;
    //     exit(1);
    // }

    // each thread will sieve equal sized subintervals
    int nthreads = 4;
    mpz_class min = sqrt_n - range, max = sqrt_n + range;
    mpz_class interval = (max-min)/nthreads;

    std::vector<std::thread> threads_vec;

    for (int i=0; i<nthreads; ++i){
        mpz_class lo = min + i*interval;
        mpz_class hi = lo + interval;
        threads_vec.emplace_back(
            smooth_sieve,
            n,lo,hi,
            chunk_size,
            std::cref(factor_base),
            std::cref(sqrts));
    }
    for (auto& t: threads_vec){
        t.join();
    }

    return 0;
}
