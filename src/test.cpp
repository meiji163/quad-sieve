#include "sieve.h"
#include "num_theory.h"
#include <iostream>

void print( const fact_map& f){
	for (auto [p,e] : f){
		std::cout << p << "^" << e << " ";
	}
	std::cout << std::endl;
}

void test_factor(){
	auto fac = factor(1234321);
	for (auto [p,e]: fac){
		std::cout << p << "^" << e << " ";
	}
	std::cout << std::endl;
}


void test_smooth_sieve(){
	int64_t n = 1234321;
	int64_t B = 1000;
	auto primes = prime_sieve(B);
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
	auto cand = smooth_sieve(2*B,10,sprm,sqrts);
	std::cout << "Candidates: " << std::endl;
	for (auto c: cand){
		std::cout << c << " ";
	}
	std::cout << std::endl;

	auto found = find_smooth(n, B, cand, sprm, primes);
	std::cout << "Smooth Numbers: " << std::endl;
	for ( const auto [p,f] : found){
		std::cout << p << ": ";
		print(f);
	}
}

int main(){
	test_smooth_sieve();
	return 0;
}
