#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <map>

#pragma once

/* n^e (mod p) */
int64_t mod_pow(int64_t n, int64_t e, int64_t p);

/* check if n is a square (quadratic residue) mod p */
bool is_residue(int64_t n, int64_t p);

/* find a square root of n mod p, or return -1 if it doesn't exist*/
int64_t mod_sqrt(int64_t n, int64_t p);

/* trial division factorization */
std::map<int64_t,int> factor(int64_t n);

/* factor with base of primes, returning vector of exponents */
std::vector<int> factor_with_base(mpz_class n, const std::vector<int64_t> &base);
