CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic

all: test_matmodel

test_matmodel: matmodel.cpp matmodel.h matrix.cpp matrix.h vector.cpp vector.h test_matmodel.cpp
	$(CXX) $(CXXFLAGS) vector.cpp matrix.cpp matmodel.cpp test_matmodel.cpp -o $@

clean:
	rm -f test_matmodel *.o
