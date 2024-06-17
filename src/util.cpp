#include "util.h"

int64_t mpz_to_int64(mpz_class z)
{
    unsigned long long result = 0;
    mpz_export(&result, 0, -1, sizeof(result), 0, 0, z.get_mpz_t());
    return (int64_t)result;
}
