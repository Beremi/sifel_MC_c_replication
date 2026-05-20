CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic

all: test_mohrc_ugn compare_matlab_export test_mohrc3d_ugn compare_matlab_export_3d

test_mohrc_ugn: mohrc_ugn.cpp mohrc_ugn.h iotools.h global.h mechmat.h sifel_compat.cpp matrix.cpp matrix.h vector.cpp vector.h test_mohrc_ugn.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp sifel_compat.cpp mohrc_ugn.cpp test_mohrc_ugn.cpp -o $@

compare_matlab_export: mohrc_ugn.cpp mohrc_ugn.h iotools.h global.h mechmat.h sifel_compat.cpp matrix.cpp matrix.h vector.cpp vector.h compare_matlab_export.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp sifel_compat.cpp mohrc_ugn.cpp compare_matlab_export.cpp -o $@

test_mohrc3d_ugn: mohrc3d_ugn.cpp mohrc3d_ugn.h iotools.h global.h mechmat.h sifel_compat.cpp matrix.cpp matrix.h vector.cpp vector.h test_mohrc3d_ugn.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp sifel_compat.cpp mohrc3d_ugn.cpp test_mohrc3d_ugn.cpp -o $@

compare_matlab_export_3d: mohrc3d_ugn.cpp mohrc3d_ugn.h iotools.h global.h mechmat.h sifel_compat.cpp matrix.cpp matrix.h vector.cpp vector.h compare_matlab_export_3d.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp sifel_compat.cpp mohrc3d_ugn.cpp compare_matlab_export_3d.cpp -o $@

clean:
	rm -f test_mohrc_ugn compare_matlab_export test_mohrc3d_ugn compare_matlab_export_3d *.o
