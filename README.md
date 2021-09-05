# Quadratic Sieve

## Outline
Input: BIG int n we want to factor

### Part 1: Sieving
1. Find factor base F: 
	primes below smoothness bound B such that n is a square mod p 

2. Sieve [ x^2 - n | x in [-M,M] ] to find F-smooth numbers.
	a^2 = n (mod p)  ==> (a + k*p)^2 = n (mod p)

3. Factorize each F-smooth number 

### Part 2: Linear Alg

1. Get mod2 exponent vectors from the F-smooth numbers 
	and find a linear dependence (over Z/Z2)

Need to find kernel of |F| x K matrix 
where K is number of F-smooth number found.

Techniques:
	- Gaussian elimination: O(K^3), K > 1e5 too slow 
	- Structured Gaussian elimination, K>1e5 too much space 
	- Block Lanczos
	- Block Wiedemann (using Berlekamp-Massey)

Finally get a factorization 
x1^2 * x2^2 * ... xm^2 = (x1^2 - n) * (x2^2-n) * ... (xm^2 - n)

where both sides are squares. So x^2 = y^2 (mod n), (x-y)(x+y) = 0 (mod n)
Hopefully this gives you a nontrivial divisor of n.
