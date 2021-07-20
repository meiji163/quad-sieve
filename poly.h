#include <set> 
#include <stdint.h>

struct Z2_Poly{
	std::set<uint64_t> coeff;
	Z2_Poly& operator+=( Z2_Poly const& p);
	Z2_Poly& operator-=( Z2_Poly const& p);
};

Z2_Poly operator+ ( Z2_Poly const& p1, Z2_Poly const& p2);
Z2_Poly operator- ( Z2_Poly const& p1, Z2_Poly const& p2);
Z2_Poly operator* ( Z2_Poly const& p1, Z2_Poly const& p2);

