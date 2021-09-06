#include <boost/multiprecision/cpp_int.hpp>
#include <set>
#include <vector>
#include <cassert>

struct SparseMat{
	uint64_t rows;
	uint64_t cols;

	std::map<uint64_t, set<uint64_t> entry_map;
	uint64_t rank;
};	

std::set<std::vector<uint64_t>> GetKernelVectors(SparseMat M);

bool dot(const std::vector<bool>& v1, const std::vector<bool>& v2);
