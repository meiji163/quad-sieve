#include "sieve.h"
#include "num_theory.h"
#include <iostream>

void test_factor(){
	auto fac = factor(1234321);
	for (auto [p,e]: fac){
		std::cout << p << "^" << e << " ";
	}
	std::cout << std::endl;
}


void test_smooth_sieve(){
	auto primes = prime_sieve(100);
	int64_t n = 1029;
	std::vector<int64_t> sqrts;
	std::vector<int64_t> sprm;
	for (int i=1; i<primes.size(); ++i){
		int64_t s = mod_sqrt(n, primes[i]);
		if (s>=0){
			sqrts.push_back(s);
			sprm.push_back(primes[i]);
			std::cout << primes[i] << "  " << s << std::endl;
		}
	}
	auto lgs = smooth_sieve(n,100,sprm,sqrts);
}

int main(){
	test_factor();
	return 0;
}
