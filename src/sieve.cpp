#include "sieve.h"
#include <iostream>

typedef unsigned long long ull;

int ilog2(int64_t n){
    if (n <=0){
        throw std::runtime_error("ilog2 domain error"); 
    }
    int out = 0;
    while (n >>= 1) ++out;
    return out; 
}

std::vector<int64_t> prime_sieve(ull max){
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

std::vector<mpz_class> smooth_sieve(
        mpz_class n,
        mpz_class min,
        mpz_class max,
        const std::vector<int64_t>& factor_base,
        const std::vector<int64_t>& sqrts)
{
    int64_t B = 100000000; // chunk size
    int toler = 30; // tolerance for searching for approximate smooth numbers

    std::vector<int64_t> offset(factor_base.size(),0);
    std::vector<int64_t> offset_neg(factor_base.size(),0);
    std::vector<int> logs(B,0);
    std::vector<int> ilog2_vec(factor_base.size(),0); // precomputed ilog2()
    for (int i=0; i<factor_base.size(); ++i){
        ilog2_vec[i] = ilog2(factor_base[i]);
    }

    mpz_class L = min;
    if (L%2==1) L++;

    // index i <--> L + 2*i + 1
    // L is even
    //

    // initialize offsets
    // L+1 + 2*offset[i] = sqrts[i] (mod p[i])
    // L+1 + 2*offset_neg[i] = p - sqrts[i] (mod p[i])
    int64_t p, s, j;
    int i;
    std::cout << "L=" << L << std::endl;
    for (int i=0; i<factor_base.size(); ++i){
        p = factor_base[i];
        s = sqrts[i];
        mpz_class a = (s-L-1)%p;
        while (a<0) a += p;

        j = mpz_to_int64(a);
        if (j%2==1) j += p;
        j /= 2;
        offset[i] = j;

        assert((L+1+2*j)%p == s);

        a = (p-s-L-1)%p;
        while (a<0) a += p;
        j = mpz_to_int64(a);
        if (j%2==1) j += p;
        j /= 2;
        offset_neg[i] = j;

        assert((L+1+2*j)%p == p-s);
    }
    std::cout << "[INFO] start smooth sieve" << std::endl;


    std::vector<mpz_class> smooth; // candidate smooth
    int log_p;
    while(L < max){
        for (int i=0; i<factor_base.size(); ++i){
            p = factor_base[i];
            log_p = ilog2_vec[i];
            for (j=offset[i]; j<B; j+=p) {
                logs[j] += log_p;
            }
            for (j=offset_neg[i]; j<B; j+=p) {
                logs[j] += log_p;
            }
        }

        for (int64_t i=0; i<logs.size(); ++i){
            mpz_class x = L+1 + 2*i;
            mpz_class y = x*x - n;
            int64_t target_log = mpz_sizeinbase(y.get_mpz_t(), 2);
            if (logs[i] + toler > target_log){
                std::cout << logs[i] << " " << target_log << " "
                          << x << std::endl;
                smooth.push_back(x);
            }
        }

        // shift to next chunk
        for (int i=0; i<offset.size(); ++i){
            p = factor_base[i];
            offset[i] = (offset[i]-B)%p;
            while(offset[i]<0) offset[i]+=p;

            offset_neg[i] = (offset_neg[i]-B)%p;
            while(offset_neg[i]<0) offset_neg[i]+=p;
        }
        for (int i=0; i<logs.size(); ++i){
            logs[i] = 0;
        }
        L = L+2*B;
    }
    return smooth;
}

bool _isin(int64_t d, const std::vector<int64_t>& V){
    auto it = std::lower_bound(V.begin(),V.end(),d);
    return (it!= V.end() && *it == d);
}

std::map<int64_t,fact_map> find_smooth(
        int64_t n, int64_t B,
        const std::vector<int64_t>& cds, 
        const std::vector<int64_t>& facb,
        const std::vector<int64_t>& prms)
{
    std::map<int64_t,fact_map> smth;

    //large prime variation
    // x^2 - n = p1 .. pk * L for L a large prime
    struct PartialFacts{ 
      // L -> [ x .. ]
      std::map<int64_t, std::vector<int64_t> > Lmap;

      // x -> partial factorization of x^2 - n 
      std::unordered_map<int64_t, fact_map> facts; 
    } partial;

    for (auto i=cds.begin(); i!=cds.end(); ++i){
        fact_map facs;
        int64_t c = (*i) * (*i) - n;

        if ( c == 0 ){
            //perfect square (should check for this earlier)
            facs[*i] = 2;
            smth[n] = facs;
            return smth;
        } else if ( c < 0 ){
            facs[-1] = 1;
            c = -c;
        }

        bool brk = false;
        for (auto j=prms.begin(); j != prms.end(); ++j){
            if ( c % (*j) == 0){
                if ( _isin(*j,facb) || *j == 2) {
                    do {
                        facs[ *j ]++;
                        c /= (*j);
                    } while ( c % (*j) == 0);
                }else{
                    // divisible by prime < B not in factor base
                    brk = true;
                    break;
                }
            }
        }

        if (c == 1){
            smth[*i] = facs;
        }else if (!brk && c > B && c < 20*B){
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

                //std::cout << x0 << " * " << xi << " = " << x0*xi << std::endl;

                //factorization of (x0^2 - n) * (xi^2 - n)
                for (auto [p,e]: partial.facts[x0]) 
                    facs[p] += e;
                facs[L] = 2;
                smth[x0*xi] = facs;
            }
        }
    }
    return smth;
}
