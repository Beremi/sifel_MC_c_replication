CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic

all: test_mcppp2d compare_matlab_export

test_mcppp2d: mcppp2d_core.cpp mcppp2d_core.h test_mcppp2d.cpp
	$(CXX) $(CXXFLAGS) mcppp2d_core.cpp test_mcppp2d.cpp -o $@

compare_matlab_export: mcppp2d_core.cpp mcppp2d_core.h compare_matlab_export.cpp
	$(CXX) $(CXXFLAGS) mcppp2d_core.cpp compare_matlab_export.cpp -o $@

clean:
	rm -f test_mcppp2d compare_matlab_export *.o
