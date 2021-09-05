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

std::vector<int64_t> smooth_sieve(int64_t M, int64_t tol, 
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

	std::vector<int64_t> cand;
	for(int64_t i=0; i<lgs.size(); ++i){
		if( lgs[i] >= ilog2( std::abs( 2*(i-M) + 1)) + tol ){
			cand.push_back(2*(i-M)+1);
		}
	}
	return cand;
}

bool _isin(int64_t d, const std::vector<int64_t>& V){
	auto it = std::lower_bound(V.begin(),V.end(),d);
	return it!= V.end();
}

std::map<int64_t,fact_map> find_smooth(
		int64_t n, int64_t B,
		const std::vector<int64_t>& cds, 
		const std::vector<int64_t>& facb,
		const std::vector<int64_t>& prms)
{
	std::map<int64_t,fact_map> smth;

	//large prime variation
	// x^2 - n = p1 .. pk * L 
	struct PartialFacts{ 
		// L --> [ x .. ]
		std::map<int64_t, std::vector<int64_t> > Lmap;

		// x --> partial factorization of x^2 - n 
		std::unordered_map<int64_t, fact_map> facts; 
	} partial;


	for (auto i=cds.begin(); i!=cds.end(); ++i){
		fact_map facs;
		int64_t c = (*i) * (*i) - n;
		for (auto j=prms.begin(); j != prms.end(); ++j){
			if ( c % (*j) == 0 && _isin(*j,facb)){ 
				do {
					facs[ *j ]++;
					c /= *j;
				} while ( c % (*j) == 0);
			}else{
				// divisible by prime < B not in factor base
				break;
			}
		}

		if (c == 1){
			smth[*i] = facs;
		}else if (c != 1 && c > B && c < 30*B){
			partial.Lmap[c].push_back(*i);
			partial.facts[*i] = facs;
		}
	}

	for (auto [L,v]: partial.Lmap){
		if (v.size() >1){
			int64_t x0 = v[0];
			for(int i=1; i<v.size(); ++i){
				int64_t xi = v[i];
				fact_map facs = partial.facts[xi];

				//factorization of x0 * xi 
				for (auto [p,e]: partial.facts[x0]) 
					facs[p] += e;
				facs[L] = 2;
				
				smth[x0*xi] = facs;
			}
		}
	}
	return smth;
}
