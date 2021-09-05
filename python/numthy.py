import math
import numpy as np
from numba import njit

@njit
def mod_pow(n, e, p):
    '''n^e mod p'''
    exp = e
    out = 1
    pw = n%p
    while exp > 0:
        if exp & 1:
            out = (out*pw) %p
        pw = (pw*pw)%p
        exp >>=1
    return out%p

@njit
def is_residue(n, p):
    if mod_pow(n, (p-1)//2, p) == 1:
        return True
    return False

@njit
def mod_sqrt(n, p):
    '''
    Find a square root of n modulo prime p using Shanks-Tonelli algorithm.
    Returns:
        sqrt in [0, p-1] or -1 if it does not exist.
    '''
    if p==2:
        return n%2
    elif p%4 == 3:
        root = mod_pow(n%p, (p+1)//4, p)
        return root if (root*root)%p == n%p else -1
    # p = 1 mod 4
    # p-1 = Q * 2^S for Q odd
    S = 0
    Q = p-1
    while Q%2 == 0:
        Q = Q//2
        S += 1
    non_res = 2
    while is_residue(non_res%p, p):
        non_res += 1
    c = mod_pow(non_res%p, Q, p)
    t = mod_pow(n%p, Q, p)
    root = mod_pow(n%p, (Q+1)//2, p)

    # c^( 2^(S-1) ) = -1 mod p
    # t^( 2^(S-1) ) = 1 mod p
    # root^2 = t*n mod p
    while t>1:
        e = 0
        pw = t
        while pw != 1 and e<= S:
            pw *= pw
            pw %= p
            e += 1
        if e>= S:
            return -1 #not a residue
        b = mod_pow(c%p, 1 << (S-e-1), p)
        S = e
        c = (b*b) % p
        t = (t*c) % p
        root = (root*b) % p
    return root if t ==1 else 0

@njit
def sieve_eratosthenes(M):
    '''Sieve for all primes less than M
    returns:
        list of primes
    '''
    is_prime = [True]*((M+2)//2)
    primes = [2,3] 
    i = 2
    p = 5
    is_2mod3 = True
    while p<M:
        if is_prime[i]:
            b = is_2mod3
            primes.append(p)
            j = 2*i*(i+1)
            while 2*j+1 < M:
                is_prime[j] = False
                incr = 1 if b else 2
                j += p*incr
                b = not b
        incr = 1 if is_2mod3 else 2
        i += incr
        p += 2*incr
        is_2mod3 = not is_2mod3
    return primes

def log_sieve(n, fbase, logs):
    roots = dict()
    sqrt_n = math.ceil(math.sqrt(n))
    M = (len(logs)-1)//2
    for i in range(len(fbase)):
        p = fbase[i]
        res = mod_sqrt(n%p, p)
        roots[p] = res
        if res < 0:
            continue
        for r in (res,p-res):
            low = (r+M-sqrt_n) % p 
            logs[np.arange(low,2*M+1,p)] -= math.ceil(math.log2(p))
    return roots

def try_factor(x, fbase):
    '''Try to factor x over the factor base
    params:
        x (int) : number to factor
        fbase (list): list of primes 

    returns: 
        dict of factors with multiplicity and,
        remainder ( =1 if x factored completely)
    '''
    facts = dict()
    facts[-1] = 1 if x < 0 else 0
    x = abs(x)
    for p in fbase:
        while x%p == 0:
            x = x//p
            if p not in facts:
                facts[p] = 0
            facts[p] += 1
    return facts, x

def div8_trans(n):
    ''' Multiply n by constant so that x^2 - n is divisible by 8 when it is even '''
    r = n%8
    if r == 3:
        return 5*n
    elif r == 5:
        return 3*n
    elif r == 7:
        return 7*n
    else:
        return n

def qs_bound(n):
    '''The asymptotic bound for the quadratic sieve'''
    logn = math.log(n)
    return math.exp(math.sqrt(logn * math.log(logn)))
