import math
import numpy as np
from scipy import sparse

def berlekamp_massey(seq):
    ''' Find the minimal polynomial for the binary sequence.
    Returns:
        set of representing polynomial; 
        n is in the set if and only if the coeff of x^n is 1
    '''
    N = len(seq)
    s = seq[:]
    k = 0
    while(s[k] == 0):
        k += 1
        if k >= N:
            raise ValueError 
    f = set([k + 1, 0])
    l = k + 1
    g = set([0])
    a = k
    b = 0
    for n in range(k+1, N):
        d = 0
        for ele in f:
            d ^= s[ele + n - l]
        if d == 0:
            b += 1
        else:
            if 2 * l > n:
                f ^= set([a - b + ele for ele in g])
                b += 1
            else:
                temp = f.copy()
                f = set([b - a + ele for ele in f]) ^ g
                l = n + 1 - l
                g = temp
                a = b
                b = n - l + 1
    return f

def make_Z2mat(smooth, fbase):
    '''Make sparse matrix whose columns are the mod 2 exponents of 
    smooth ints (factorized over the factor base)

    params:
        smooth: list of dicts containing the factorizations
        fbase: the factor base

    returns: 
        factor matrix (scipy.sparse.dok_matrix)
    '''
    fb_idx = {p:i for i, p in enumerate(fbase)}
    fb_idx[-1] = len(fb_idx)
    mat = sparse.dok_matrix(
                (len(fbase)+1, len(smooth)),
                dtype = np.uint64
            )
    for i in range(len(smooth)):
        for p in smooth[i]:
            mat[fb_idx[p], i] = 1
    return mat

def sym_square(mat, v):
    ''' (mat.T @ mat) . v mod 2 ''' 
    return mat.transpose().dot( mat.dot(v)%2 )%2

def linear_seq(mat, u, v, seq_len):
    ''' Compute u.T . (mat.T @ mat)^k . v for k=0,...,seq_len.

    params:
        mat: N x M scipy sparse matrix
        u, v : length M  numpy arrays 
        seq_len (int): the length of the sequence

    returns:
        binary sequence (list)
    '''
    seq = []
    w = np.copy(v)
    for _ in range(seq_len):
        seq.append( bool(u.dot(w)%2) )
        w = sym_square(mat, w)
    return seq

def eval_poly(poly, mat, v):
    ''' Plug the square mat.T @ mat into polynomial and evaluate it on a vector.

    params:
        poly (set): the mod2 poly (represented by set of nonzero coeffs)
        mat : N x M scipy sparse matrix
        v : length M numpy array
    
    returns:
        length N numpy array
    '''
    out = np.zeros_like(v)
    w = np.copy(v)
    for i in range(max(poly)+1):
        if i in poly:
            out = (out+w)%2
        w = sym_square(mat, w)
    return out%2

def weidemann(mat, seq_len = None):
    ''' Attempt to find a null vector for matrix using Wiedemann algorithm.
    
    params:
        mat: N x M scipy sparse matrix
        seq_len: length of sequence to use in Berlekamp-Massey
                (default: 2*M + 3)
    returns:
        candidate null vector for mat
    '''
    u = np.random.choice( a=[0,1], p=[0.1,0.9], size=mat.shape[1])
    v = np.random.choice( a=[0,1], p=[0.2,0.8], size=mat.shape[1]) 

    if seq_len is None:
        seq_len = 2*mat.shape[1] + 5
    seq = linear_seq(mat, u, v, seq_len) 
    poly = berlekamp_massey(seq)
    delta = min(poly)
    red_poly = {i-delta for i in poly}
    null_v = eval_poly(red_poly, mat, v)

    for i in range(delta):
        w = sym_square(mat, null_v)
        if not np.any(w):
            return null_v 
        null_v = w 
    return null_v

