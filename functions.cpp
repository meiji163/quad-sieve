#include <stdint.h>

/* calculate n^e mod p */
uint64_t mod_pow(uint64_t n, uint64_t e, uint64_t p){
	uint64_t exp = e;
	uint64_t out = 1; 
	uint64_t pow = n;
	while(exp > 0){
		if (exp & 1){
			out = (out*pow) % p;
		}
		pow = (pow*pow) % p;
		exp >>= 1;
	}
	return out;
}

/* check if n is a quadratic residue mod p */
bool is_residue(uint64_t n, uint64_t p){
	if (mod_pow(n, (p-1)/2, p)%p == 1){
		return true;
	}else{
		return false;
	}
}

/* find sqrt of n mod p */
uint64_t mod_sqrt(const uint64_t n, const uint64_t p){
	if ( p == 2){
		return (n%2);
	}
	else if ( p%4 == 3){
		uint64_t root = mod_pow(n, (p+1)/4, p);
		if ( (root*root)%p == n%p ){
			return root;
		}else{
			return -1;
		}
	}else{ /* p = 1 mod 4 */
		/* p-1 = Q*2^S */
		uint64_t S = 0;
		uint64_t Q = p-1;
		while(Q%2 == 0){
			Q /= 2;
			++S; 
		}
		uint64_t non_res = 2;
		while(is_residue(non_res, p)){
			++non_res;
		}
		uint64_t b, e, pow;
		uint64_t c = mod_pow(non_res, Q, p);
		uint64_t t = mod_pow(n, Q, p);
		uint64_t root = mod_pow(n, (Q+1)/2, p);

		/* c^( 2^(max-1) ) = -1 mod p
		 * t^( 2^(max-1) ) = 1 mod p
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
}
