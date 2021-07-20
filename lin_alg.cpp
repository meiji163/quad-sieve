#include "lin_alg.h"

std::vector<bool> operator* (const Sparse_Z2_Mat& M, const std::vector<bool>& v){
	std::assert ( v.size() == M.cols);
	std::vector<bool> out(M.rows, 0);
	std::set<uint64_t>::iterator it;
	bool b = 0;
	uint64_t r = 0;
	for (it = M.nzero.begin(); it != M.nzero.end(); ++it){
		while( *it >= r * M.rows){
			out[r++] = b; 
			b = 0;
		}
		if ( v[ (*it)% M.rows] ){
			b = !b;
		}
	}
	return out;
}

bool dot(const std::vector<bool>& v1, const std::vector<bool>& v2){
	std::assert ( v1.size() == v2.size());
	bool out = 0;
	for (int i = 0; i< v1.size(); ++i){
		out ^= v1[i] & v2[i];
	}
	return out;
}
