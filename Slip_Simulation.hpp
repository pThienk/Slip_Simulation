#ifndef SLIP_SIMULATION_HPP
#define SLIP_SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdint>
#include "cxxopts.hpp"

// Enum for print options
enum class PRINT_TYPE {
	BOTH = 0,
	STRESS_ONLY = 1,
	STRAIN_ONLY = 2
};

inline void print_to_file();

#ifdef _OLD_CODE

double* ranmarin(int ijkl, int N);

double* ranwbl(int ijkl, int N, double k, double lambda);

#endif // _OLD_CODE


#endif // !SLIP_SIMULATION_HPP
