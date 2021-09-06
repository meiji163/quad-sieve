#include <boost/multiprecision/cpp_int.hpp>
#include <set>
#include <vector>
#include <cassert>

struct SparseMat{
	uint64_t rows;
	uint64_t cols;

	std::map<uint64_t, set<uint64_t> entry_map;
	std::map<uint64_t, set<uint64_t> reverse_map;
	uint64_t rank;
};

struct F2Vector{
	uint64_t size;
	std::set<uint64_t> entries;
}

std::set<std::vector<uint64_t>> GetKernelVectors(SparseMat M);

F2Vector operator*(const SparseMat&, const F2Vector&);
