#include "sieve.h"
#include "num_theory.h"
#include <gtest/gtest.h>

const int64_t PR100[] = {2,3,5,7,11,13,17,
                        19,23,29,31,37,41,
                        43,47,53,59,61,67,
                        71,73,79,83,89,97};

TEST( PowTest, Fermat){
  for (int i=1; i<25; ++i){
    int64_t p = PR100[i];
    for (int64_t n=0; n<p; ++n){
      EXPECT_EQ( mod_pow(n, p, p), n);
    }
  }
}

TEST(SqrtTest, 1mod4){
  int64_t m = mod_sqrt( 10000, 113 );
  EXPECT_EQ( (m*m) % 113, 10000 % 113);

  m = mod_sqrt(10000, 613);
  EXPECT_EQ ( (m*m) % 613, 10000 % 613);

  m = mod_sqrt(10001, 613);
  EXPECT_EQ ( m, -1);

  m = mod_sqrt( 7*613, 613);
  EXPECT_EQ( m, 0);
}

TEST(SqrtTest, 3mod4){
  int64_t m = mod_sqrt(10000, 607);
  EXPECT_EQ( (m*m) % 607, 10000 % 607);

  m = mod_sqrt(10000, 911);
  EXPECT_EQ( (m*m) % 911, 10000 % 911);

  m = mod_sqrt(10001, 911);
  EXPECT_EQ( m, -1);
} 

TEST(FactorTest, BasicAssertions){
  auto fac = factor(1234321);
  EXPECT_EQ(fac.size(),2);
  auto it = fac.begin(); 
  EXPECT_EQ(it->first, 11);
  EXPECT_EQ(it->second, 2);
  ++it;
  EXPECT_EQ(it->first, 101);
  EXPECT_EQ(it->second, 2);
}

TEST(FactorTest, NonPositive){
	auto fac = factor(-1234321);
  EXPECT_EQ(fac.size(),3);

	auto it = fac.begin(); 
  EXPECT_EQ(it->first, -1);
  EXPECT_EQ(it->second, 1);
  ++it;
  EXPECT_EQ(it->first, 11);
  EXPECT_EQ(it->second, 2);
  ++it;
  EXPECT_EQ(it->first, 101);
  EXPECT_EQ(it->second, 2);

  fac = factor(0);
  EXPECT_TRUE(fac.empty());
}

TEST(EratostheneTest, BasicAssertion){
  auto ps = prime_sieve(100);
  for (int i=0; i<ps.size(); ++i)
    EXPECT_EQ( ps[i], PR100[i]);
}

TEST(SmoothTest, BasicAssertions){
	int64_t n = 1013*997;
	int64_t B = 100;
  int64_t M = 2*B;

	auto prm = prime_sieve(B);
	std::vector<int64_t> sqrts;
	std::vector<int64_t> sprm;

	for (int i=1; i<prm.size(); ++i){
		int64_t s = mod_sqrt(n, prm[i]);
		if (s>=0){
			sqrts.push_back(s);
			sprm.push_back(prm[i]);
		}
	}

  int64_t sqrt_n = (int64_t)std::sqrt(n);
  if (sqrt_n%2 ==1){
      sqrt_n++;
  }

    // find smooth numbers
	auto cand = smooth_sieve(n,M,10,sprm,sqrts);
	auto found = find_smooth(n, B, cand, sprm, prm);

  for( auto [x,f] : found){
    if ( x > sqrt_n + 2*M)
      break;
    int64_t c = 1;
    for(const auto [p,e] : f){
      for (int i=0; i<e; ++i){
          c *= p;
      }
    }
    EXPECT_EQ( x*x - n, c);
  }
}
