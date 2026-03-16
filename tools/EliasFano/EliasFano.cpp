#include <string>
#include <fstream>
#include <iostream>
//#include <format>
#include "2i_bench/include/util.hpp"
#include "2i_bench/include/ds2i/strict_elias_fano.hpp"
#include "2i_bench/external/succinct/bit_vector.hpp"

void int64ToChar(char a[], int64_t n) {
  memcpy(a, &n, 8);
}

int main(int argc, char *argv[]) {
	if (argc != 1) {
		std::cerr << "Compress the given mask using Elias-Fano encoding" << std::endl;
		std::cerr << "Outputs the Elias-Fano compressed representation." << std::endl;
		std::cerr << "The input is expected to be a single string of 0 a 1 on the standard input";
		return std::string(argv[1]) == "-h";
	}

	std::string input; std::cin >> input;
	std::vector<uint64_t> vec;
    vec.push_back(0);
    bool last = input[0] == '1';
	for (size_t i = 1; i < input.size(); ++i) {
        bool one = input[i] == '1';
        if (last != one) {
            vec.push_back(i);
        }
        last = one;
    }
    //std::cout << input.size() << std::endl;
    //std::cout << vec.size() << std::endl;
    //for (const auto &elem : vec) {
    //    std::cout << elem << " ";
    //}
    //std::cout << std::endl;
    succinct::bit_vector_builder bvb;
    ds2i::global_parameters params;
    ds2i::strict_elias_fano::write(bvb, vec.begin(), input.size(), vec.size(), params);
    succinct::bit_vector bv(&bvb);
    //typename strict_elias_fano::enumerator r(bv, 0, universe, vec.size(), params);
    for (size_t i = 0; i < bv.data().size(); i++) {
        uint64_t x = bv.data()[i];
        for (int j = 0; j < 8; ++j) {
            std::cout << (char)(x & ((1 << 8) - 1));
            x >>= 8;
        }
        //std::cout << std::format("{:b}", bv[i]);
    }
    //std::cout << reinterpret_cast<const char*>(bv.data().data());
    //std::cout << std::endl;
    //for (const auto &elem : bv.data()) {
    //   std::cout << elem << ",";
    //}
    //std::cout << std::endl;
	return 0;
}
