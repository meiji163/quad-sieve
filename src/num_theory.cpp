#include "num_theory.h"

int64_t mod_pow(int64_t n, int64_t e, int64_t p){
	int64_t exp = e;
	int64_t out = 1; 
	int64_t pow = n;
	while(exp > 0){
		if (exp & 1){
			out = (out*pow) % p;
		}
		pow = (pow*pow) % p;
		exp >>= 1;
	}
	return out;
}

bool is_residue(int64_t n, int64_t p){
	return (mod_pow(n, (p-1)/2, p)%p == 1);
}

// Tonelli-Shanks algorithm
int64_t mod_sqrt(int64_t n, int64_t p){
	if (p == 2){
		return (n%2);
	}
	if (p%4 == 3){
		int64_t root = mod_pow(n, (p+1)/4, p);
		if ( (root*root)%p == n%p ){
			return root;
		}else{
			return -1;
		}
	}

	/* p = 1 mod 4 */
	/* factor p-1 = Q*2^S */
	int64_t S = 0;
	int64_t Q = p-1;
	while(Q%2 == 0){
		Q /= 2;
		++S;
	}

	int64_t non_res = 2;
	while(is_residue(non_res, p)){
		++non_res;
	}
	int64_t b, e, pow;
	int64_t c = mod_pow(non_res, Q, p);
	int64_t t = mod_pow(n, Q, p);
	int64_t root = mod_pow(n, (Q+1)/2, p);

	/* The loop maintains the invariants:
	 * c^{ 2^{max-1} } = -1 mod p
	 * t^{ 2^{max-1} } = 1 mod p
	 * root^2 = t*n mod p */
	while(t>1){
		e = 0;
		pow = t;
		while( pow != 1 && e<= S){
			pow *= pow;
			pow %= p;
			++e;
		}
		if ( e >= S){
			return -1;
		}
		b = mod_pow(c, 1 << (S-e-1), p);
		S = e;
		c = (b*b) % p;
		t = (t*c) % p;
		root = (root*b) % p;
	}
	if (t == 1){
		return root;
	}else{
		return 0;
	}
}

std::map<int64_t,int> factor(int64_t n){
	std::map<int64_t,int> facs;
    if ( n == 0){
        return facs;
    } else if (n < 0){
        facs[-1] = 1;
        n = -n;
    }

	while(n%2 == 0){
		facs[2]++;
		n /= 2;
	}
	for (int64_t d=3; d*d <= n; d+=2){
	  while( n%d == 0){
			facs[d]++;
			n /= d;
		}
	}
	if (n != 1) facs[n]++;
	return facs;
}

std::vector<int> factor_with_base(mpz_class n, const std::vector<int64_t> &base){
	std::vector<int> exp(base.size(), 0);

	for (int i=0; i<base.size(); ++i){
		int64_t p = base[i];
		int e = 0;
		while(n%p==0){
			e++;
			n = n/p;
		}
		exp[i] = e;
	}
	return exp;
}
