#ifndef UTIL_H_
#define UTIL_H_

#include "gmp.h"
#include "gmpxx.h"

/* convert mpz_class to int64_t */
int64_t mpz_to_int64(mpz_class z);

#endif // UTIL_H_
