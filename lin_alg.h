#include <boost/multiprecision/cpp_int.hpp>
#include <set>
#include <vector>
#include <cassert>

struct Sparse_Z2_Mat{
	uint64_t rows;
	uint64_t cols;

	/* store the nonzero coordinates (i,j) as rows*i + j */
	std::set<uint64_t> nzero;
};

std::vector<bool> operator* (const Sparse_Z2_Mat& M, const std::vector<bool>& v);

bool dot(const std::vector<bool>& v1, const std::vector<bool>& v2);
