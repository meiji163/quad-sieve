#include <bitset>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdint.h>

#pragma once 
typedef std::map<int64_t,int> fact_map;
int ilog2( int64_t n);
std::vector<int64_t> prime_sieve(int64_t max);
std::vector<int64_t> smooth_sieve(
			int64_t n, int64_t M, int64_t tol,
			const std::vector<int64_t>& facb, 
			const std::vector<int64_t>& sqrts
		);

std::map<int64_t,fact_map> find_smooth(int64_t n, int64_t B,
			const std::vector<int64_t>& cds,
			const std::vector<int64_t>& facb,
			const std::vector<int64_t>& prms
		);

