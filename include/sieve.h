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
#include "util.h"

#pragma once 
typedef std::map<int64_t,int> fact_map;
uint16_t ilog2(int64_t n);
std::vector<int64_t> prime_sieve(unsigned long long max);
void smooth_sieve(
    mpz_class n,
    mpz_class min,
    mpz_class max,
    uint64_t chunk_size,
    const std::vector<int64_t>& factor_base,
    const std::vector<int64_t>& sqrts
);

std::map<int64_t,fact_map> find_smooth(
    int64_t n, int64_t B,
    const std::vector<int64_t>& cds,
    const std::vector<int64_t>& facb,
    const std::vector<int64_t>& prms
);
