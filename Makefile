CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic

all: test_mohrc_ugn compare_matlab_export

test_mohrc_ugn: mohrc_ugn.cpp mohrc_ugn.h iotools.h global.h mechmat.h sifel_compat.cpp matrix.cpp matrix.h vector.cpp vector.h test_mohrc_ugn.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp sifel_compat.cpp mohrc_ugn.cpp test_mohrc_ugn.cpp -o $@

compare_matlab_export: mohrc_ugn.cpp mohrc_ugn.h iotools.h global.h mechmat.h sifel_compat.cpp matrix.cpp matrix.h vector.cpp vector.h compare_matlab_export.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp sifel_compat.cpp mohrc_ugn.cpp compare_matlab_export.cpp -o $@

clean:
	rm -f test_mohrc_ugn compare_matlab_export *.o
