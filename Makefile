CXX = g++
SOURCES = Slip_Simulation.cpp
HEADERS = Slip_Simulation.hpp cxxopts.hpp
CXXFLAGS = -std=c++20 -Wall -O2
CXXFLAGS_D = -std=c++20 -Wall -Wextra -Werror -O0 -g -pedantic



slip_sim_sandbox: $(SOURCES) $(HEADERS)
	$(CXX) $(SOURCES) -o $@ -I ./ $(CXXFLAGS)
	
slip_sim_cluster: $(SOURCES) $(HEADERS)
	$(CXX) -DCLUSTER_BUILD $(SOURCES) -o $@ -I ./ $(CXXFLAGS)
	
slip_sim_debug: $(SOURCES) $(HEADERS)
	$(CXX) -DCLUSTER_BUILD $(SOURCES) -o $@ -I ./ $(CXXFLAGS_D)
	
	
	
.PHONY: clean

clean:
	rm -f slip_sim*
