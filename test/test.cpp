#include "sieve.h"
#include "num_theory.h"
#include <iostream>

void print( const fact_map& f){
	for (auto [p,e] : f){
		std::cout << p << "^" << e << " ";
	}
	std::cout << std::endl;
}

void print( const std::vector<int64_t>& v){
	for (auto n: v){
		std::cout << n << " ";
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
	int64_t n = 13923;
	int64_t B = 100;
	auto prm= prime_sieve(B);
	std::vector<int64_t> sqrts;
	std::vector<int64_t> sprm;
	for (int i=1; i<prm.size(); ++i){
		int64_t s = mod_sqrt(n, prm[i]);
		if (s>=0){
			sqrts.push_back(s);
			sprm.push_back(prm[i]);
		}
	}
	print(prm);
	print(sprm);

	auto cand = smooth_sieve(n,4*B,10,sprm,sqrts);
	std::cout << "Candidates: " << std::endl;
	for (auto c: cand){
		std::cout << c << " ";
	}
	std::cout << std::endl;

	auto found = find_smooth(n, B, cand, sprm, prm);
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
