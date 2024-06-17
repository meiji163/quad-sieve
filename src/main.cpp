#include "sieve.h"
#include "num_theory.h"

#include<gmp.h>
#include<gmpxx.h>
#include<stdint.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cassert>

typedef unsigned long long ull;

ull mpz_to_int64(mpz_class z)
{
    ull result = 0;
    mpz_export(&result, 0, -1, sizeof(result), 0, 0, z.get_mpz_t());
    return result;
}

int main(int argc, char *argv[]) {
    mpz_class n, sqrt_n;
    ull B = 1000000;
    n = "114381625757888867669235779976146612010218296721242362562561842935706935245733897830597123563958705058989075147599290026879543541";
    sqrt_n = sqrt(n);

    // std::cout << sqrt_n << std::endl
    //           << mpz_sizeinbase(sqrt_n.get_mpz_t(), 2) << std::endl;

    std::vector<int64_t> primes = prime_sieve(B);
    std::vector<int64_t> factor_base;
    std::vector<int64_t> sqrts;

    for(const auto p: primes){
        int64_t n_modp = mpz_to_int64(n%p);
        if(p>=5 && is_residue(n_modp, p)){
            int64_t root = mod_sqrt(n_modp, p);
            assert(root > 0);
            factor_base.push_back(p);
            sqrts.push_back(root);
        }
    }
    std::cout << "[INFO] " << factor_base.size() << " primes found" << std::endl;

    mpz_class min = sqrt_n - 100000000, max = sqrt_n + 100000000;
    auto smooth = smooth_sieve(n, min, max, factor_base, sqrts);
    std::cout << "[INFO] " << factor_base.size() << " candidate smooth numbers found" << std::endl;

    for (auto x: smooth) {
        mpz_class y = x*x - n;
        std::vector<int> exp = factor_with_base(y, factor_base);
        mpz_class prod = 1;

        for (int i=0; i<factor_base.size(); ++i) {
            int64_t p = factor_base[i];
            int e = p;
            while(--e){
                prod = p * prod;
            }
        }
        std::cout << prod << " " << y << std::endl;
    }

    return 0;
}
