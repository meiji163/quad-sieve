#include <bitset>
#include <vector>
#include <stdint.h>

int ilog2( int64_t n);
std::vector<int64_t> prime_sieve(int64_t max);
std::vector<int> smooth_sieve(
			uint64_t n, int64_t M, 
			const std::vector<int64_t>& facb, 
			const std::vector<int64_t>& sqrts
		);
