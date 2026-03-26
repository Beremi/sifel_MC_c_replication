CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic

all: test_matmodel compare_matlab_export

test_matmodel: matmodel.cpp matmodel.h matrix.cpp matrix.h vector.cpp vector.h test_matmodel.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp matmodel.cpp test_matmodel.cpp -o $@

compare_matlab_export: matmodel.cpp matmodel.h matrix.cpp matrix.h vector.cpp vector.h compare_matlab_export.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp matmodel.cpp compare_matlab_export.cpp -o $@

clean:
	rm -f test_matmodel compare_matlab_export *.o
