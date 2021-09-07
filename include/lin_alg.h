#include <boost/multiprecision/cpp_int.hpp>
#include <set>
#include <vector>
#include <cassert>

struct DenseMat{
	uint64_t rows;
	uint64_t cols;

	std::vector<std::vector<bool>> entries;
	DenseMat transpose(){
		DenseMat result;
		result.rows = cols;
		result.cols = rows;
		std::vector<std::vector<bool>> new_entries;
		for (uint64_t i = 0; i < cols; i++){
			std::vector<bool> new_row;
			for (uint64_t j = 0; j < rows; j++){
				new_row.push_back(entries[j][i]);
			}
			new_entries.push_back(new_row);
		}
	result.entries = new_entries;
	return result;
	}
}

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
SparseMat operator*(const SparseMat&, const SparseMat&);

