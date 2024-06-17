#include <map>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <cassert>

#include <gmp.h>
#include <gmpxx.h>
#include <stdint.h>

#pragma once 
typedef std::map<int64_t,int> fact_map;
int ilog2( int64_t n);
std::vector<int64_t> prime_sieve(unsigned long long max);
std::vector<mpz_class> smooth_sieve(
    mpz_class n,
    mpz_class min,
    mpz_class max,
    const std::vector<int64_t>& factor_base,
    const std::vector<int64_t>& sqrts
);

std::map<int64_t,fact_map> find_smooth(
			int64_t n, int64_t B,
			const std::vector<int64_t>& cds,
			const std::vector<int64_t>& facb,
			const std::vector<int64_t>& prms
		);





