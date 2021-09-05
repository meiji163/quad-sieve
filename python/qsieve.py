import math
from collections import defaultdict 
from numthy import *
from Z2alg import *

def Q(x,n):
    return (x+math.ceil(math.sqrt(n)))**2 - n

def quad_sieve(n, primes, **kwargs):
    n_t = div8_trans(n) #ensure Q(x) is divisible by 8 if it's even
    sqrt_n = math.ceil(math.sqrt(n_t))
    M = kwargs.get("M", sqrt_n//30)
    bound = kwargs.get("B", 
                    math.ceil(math.sqrt(qs_bound(n_t)))
                    )
    t = kwargs.get("t", 1.5)

    fbase = list(
                filter(lambda p: is_residue(n_t%p, p), primes) 
                )
    if len(fbase) < bound:
        raise ValueError(f"{bound} primes required for the factor base, \
                            but {len(fbase)} found")
    fbase = fbase[:bound]
    sqrs = list( 
                map( lambda x: Q(x, n_t), range(-M,M+1))
                )
    logs = list(
                map( lambda x: math.ceil(math.log2(abs(x))), sqrs)
                )
    logs = np.array(logs, dtype=np.int32)
    roots = log_sieve(n, fbase, logs)
    thresh = math.ceil(
                0.5*math.log2(n_t) + math.log2(M) - t*math.log2( fbase[-1])
                )

    # find the candidates that factor over the factor base
    cand = np.arange(-M, M+1)
    cand = cand[logs<thresh]
    np.random.shuffle(cand)
    xs = []
    smooth = []
    count = 0
    for x in cand:
        if count >= len(fbase)+2:
            break
        factors, rem = try_factor(sqrs[x+M], fbase)
        if rem == 1:
            smooth.append(factors)
            xs.append(x + sqrt_n)
            count += 1

    print(f"Found {len(smooth)} smooth numbers")

    # form matrix from mod 2 whose columns are 
    # the mod 2 exponents of factors (including -1 to account for sign)
    mat = make_Z2mat(smooth, fbase)
    mat = sparse.csc_matrix(mat)
    max_tries = 300 
    tries = 0
    null_v = np.ones(mat.shape[1])

    #try to find a null vector for mat
    while True: 
        if tries == max_tries:
            print(f"{max_tries} tries exceeded")
            return
        if np.any(mat.dot(null_v)%2):
            null_v = weidemann(mat)
        else: #success
            break
        tries += 1

    # (A+B)(A-B) = 0 mod n
    A, B = 1, 1
    B_factors = defaultdict(int)
    for i in range(len(null_v)):
        if null_v[i]:
            A = ( A * xs[i])% n_t
            for d in smooth[i]:
                B_factors[d] += smooth[i][d]

    for d, e in B_factors.items():
        if e%2 != 0:
            print(f"Bad factor: {d}")
            print(B_factors)
            return # HOW DOES THIS HAPPEN?
        if d != -1:
            B = (B * pow(d, e//2, n_t)) % n_t

    d1 = math.gcd( A-B, n_t)
    d2 = math.gcd( A+B, n_t)
    return d1, d2
