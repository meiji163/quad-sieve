#include "sieve.h"
#include <iostream>
int ilog2( int64_t n){
	if (n <=0){
		throw std::runtime_error("ilog2 domain error"); 
	}
	int16_t out = 0;
	while (n >>= 1) ++out;
	return out; 
}

std::vector<int64_t> prime_sieve(int64_t max){
	std::vector<int64_t> primes{ 2, 3 };
	std::vector<bool> is_prime(max+2,1);
	bool a, b = true;
	int64_t j, i = 5;
	while( i < max){
		if(is_prime[i]){
			primes.push_back(i);
			a = b;
			j = i*i;
			while(j< max){
				is_prime[j] = false;
				j += i*( a ? 2: 4);
				a = !a;
			}
		}
		i += b ? 2: 4;
		b = !b;
	}
	return primes;
}

std::vector<int> smooth_sieve(uint64_t n, int64_t M, 
		const std::vector<int64_t>& facb, const std::vector<int64_t>& sqrts)
{
	if (facb.size() != sqrts.size())
		throw std::runtime_error("Number of primes not equal to number of square roots");

	// sieve needs 4*M bytes 
	std::vector<int> lgs(2*M,0);

	// k <--> lgs[ (k-1)/2 + M ]
	for (int j=0; j<facb.size(); ++j){
		int64_t p = facb[j], sq = sqrts[j];
		int64_t sqs[2] = { sq, p - sq };
		for( int b=0; b<2; ++b){
			int64_t i = (sqs[b]-1)/2 + M;
			for (int64_t k=i; k<2*M; k += p){ 
				lgs[k] += ilog2(p);
			}
			for (int64_t k=i; k>=0; k-=p){
				lgs[k] += ilog2(p);
			}
		}
	}
	return lgs; 
}
