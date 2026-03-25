CXX ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra -pedantic

all: test_mcppp2d

test_mcppp2d: mcppp2d_core.cpp mcppp2d_core.h test_mcppp2d.cpp
	$(CXX) $(CXXFLAGS) mcppp2d_core.cpp test_mcppp2d.cpp -o $@

clean:
	rm -f test_mcppp2d *.o
