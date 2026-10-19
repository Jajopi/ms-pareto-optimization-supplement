.PHONY: all clean test cpptest converttest verify quick-verify pareto-test pareto-test-sanitize

CXX=         g++
# OpenMP for pareto optimization threads (disable by OMPFLAGS=), popcount instruction on x86_64
OMPFLAGS?=   -fopenmp
ifeq ($(shell uname -m),x86_64)
ARCHFLAGS?=  -mpopcnt
endif
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2 -I/opt/homebrew/include $(OMPFLAGS) $(ARCHFLAGS)
LDFLAGS=     -lz -L/opt/homebrew/lib $(OMPFLAGS)
LARGEFLAGS=  -DLARGE_KMERS
ELARGEFLAGS=  -DEXTRA_LARGE_KMERS
SRC=         src
TESTS=       tests
HEADERS=     $(wildcard $(SRC)/*.h $(SRC)/*/*.h $(SRC)/*.hpp)
GTEST=       $(TESTS)/googletest/googletest
DATA=        data


all: kmercamel

test: cpptest verify

mask-verify:
	./verify.py --k 13 --superstring_path $(DATA)/global-k13c.fa $(DATA)/spneumoniae.fa
	./verify.py --k 31 --superstring_path $(DATA)/global-k31c.fa $(DATA)/spneumoniae.fa
	./verify.py --k 63 --superstring_path $(DATA)/global-k63c.fa $(DATA)/spneumoniae.fa
	./verify.py --k 127 --superstring_path $(DATA)/global-k127c.fa $(DATA)/spneumoniae.fa


verify: verify.py kmercamel
	./verify.py --simplitigs --k 31 $(DATA)/simplitigs-k31.fa
	./verify.py $(DATA)/spneumoniae.fa

pareto-test: kmercamel
	bash $(TESTS)/pareto/run_tests.sh

pareto-test-sanitize: kmercamel-sanitize
	BIN=./kmercamel-sanitize bash $(TESTS)/pareto/run_tests.sh

kmercamel-sanitize: $(SRC)/main.cpp $(HEADERS) src/version.h
	$(CXX) $(CXXFLAGS) -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=address,undefined $(SRC)/main.cpp -o $@ $(LDFLAGS)

quick-verify: verify.py kmercamel
	./verify.py --simplitigs --k 31 $(DATA)/simplitigs-k31.fa
	./verify.py --quick $(DATA)/spneumoniae.fa

cpptest: kmercameltest kmercameltest-large kmercameltest-extra-large
	./kmercameltest
	./kmercameltest-large
	./kmercameltest-extra-large

kmercamel: $(SRC)/main.cpp $(HEADERS) src/version.h
	./create-version.sh
	$(CXX) $(CXXFLAGS) $(SRC)/main.cpp -o $@ $(LDFLAGS)
	cp kmercamel  🐫 || true

kmercameltest: $(TESTS)/unittest.cpp gtest-all.o $(HEADERS) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS)

kmercameltest-large: $(TESTS)/unittest.cpp gtest-all.o $(HEADERS) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS) $(LARGEFLAGS)

kmercameltest-extra-large: $(TESTS)/unittest.cpp gtest-all.o $(HEADERS) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS) $(ELARGEFLAGS)

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

src/version.h: src/version
	./create-version.sh


clean:
	rm -f kmercamel kmercamel-sanitize
	rm -f 🐫 || true
	rm -f kmercameltest
	rm -f kmercameltest-large
	rm -f kmercameltest-extra-large
	rm -r -f ./bin
	rm -f gtest-all.o
	rm -f src/version.h
