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
	int64_t n = 1013*997;
	int64_t B = 100;
    int64_t M = 2*B;
    std::cout << "n = " << n << std::endl;
    std::cout << "M = " << M << std::endl;

	auto prm= prime_sieve(B);
	std::vector<int64_t> sqrts;
	std::vector<int64_t> sprm;

    std::cout << "====== Factor Base ======" << std::endl;
	for (int i=1; i<prm.size(); ++i){
		int64_t s = mod_sqrt(n, prm[i]);
		if (s>=0){
			sqrts.push_back(s);
			sprm.push_back(prm[i]);
		}
	}
    std::cout << sprm.size() << " Primes" << std::endl;
	print(sprm);

	auto cand = smooth_sieve(n,M,10,sprm,sqrts);
	std::cout << "====== Candidates ======" << std::endl;
    std::cout << cand.size() << " Candidates Found" << std::endl;
	for (auto c: cand){
		std::cout << c << " ";
	}
	std::cout << std::endl;

	auto found = find_smooth(n, B, cand, sprm, prm);
	std::cout << "====== Smooth Numbers ======" << std::endl;
    std::cout << found.size()  << " Smooth Numbers Found" << std::endl;
	for ( const auto [x,f] : found){
		std::cout << "x = " << x 
            << ", x^2-n = " << x*x-n <<" = ";
		print(f);
	}
}

int main(){
	test_smooth_sieve();
	return 0;
}
