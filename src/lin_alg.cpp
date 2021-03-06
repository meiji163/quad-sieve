#include "lin_alg.h"

bool dot(const F2Vector v1, const F2Vector v2){
	std::assert (v1.size == v2.size);
	std::set<uint64_t> vec1 = v1.entries;
	std::set<uint64_t> vec2 = v2.entries;

	bool out = 0;
	for (uint64_t i: vec1){
		if (vec2.find(i) != vec2.end()) out = 1 - out;
	}
	return out;
}

F2Vector operator*(const SparseMat& M, const F2Vector& v){
	F2Vector w;
	w.size = M.rows;
	w.entries = {};
	for (int i = 0; i < w.size(); i++) {
		F2Vector row;
		row.size = M.cols;
		row.entries = M.entry_map()[i];
		if (dot(row, v)) w.insert(i);
	}
	return v;
}

SparseMat operator*(const SparseMat& M1, const SparseMat& M2){
	SparseMat R;
	//TODO: include a check for size compatability
	uint64_t a, b, c, d;
	a = M1.rows;
	b = M1.cols;
	c = M2.rows;
	d = M2.cols;
	R.rows = a;
	R.cols = d;
	for (auto {i, s} : entry_map) {
	
	}
	return R;
}
