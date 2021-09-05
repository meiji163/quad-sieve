#include "poly.h"

Z2_Poly& Z2_Poly::operator+=( Z2_Poly const& p){
	std::set<uint64_t>::iterator it, search;
	for ( it = this->coeff.begin(); it != this->coeff.end(); ++it){
		search = p->coeff.find(*it);
		if ( search != p->coeff.end()){
			this->coeff.erase(*it);
		}else{
			this->coeff.insert(*it);
		}
	}
	return *this;
}

uint64_t Z2_Poly::deg(){
	if (!this->coeff.empty()){
		return this->coeff.rbegin();
	}
	else{
		return 0;
	}
}

Z2_Poly operator+ (Z2_poly p1, Z2_Poly poly& p2){
	p1 += p2;
	return p1;
}

Z2_Poly& Z2_Poly::operator-=( Z2_Poly const& p){
	return *this += p;
}

Z2_Poly gcd( Z2_Poly const& p1, Z2_Poly const& p2){

}
