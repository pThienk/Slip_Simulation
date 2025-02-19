/*
* Remake Version: Punnatorn Thienkingkeaw (Porpun) Feb 10 2025
* [FOR DETAILS OF THE IMPROVED FEATURES, CONSULT README]
* 
* Original Created: Alan Long Jun 20 2019
*
* In service of the Dahmen Research Group's Avalanche Initiative. This code simulates our slip mean-field model exactly.
* 
* IMPORTANT!!!
*
* The system has six variables:
* (1) time_max (int) is the length of time you wish to simulate (in steps simulation timesteps).
* (2) Area (int) is the system size (in cells). Default is 1000.
* (3) consv (double) is the conservation parameter c.
* (4) rate (double) is the strain rate, set to 0 for an adiabatic system.
* (5) d (double) is the spread of arrest stresses with failure stress normalized to 1.
* (6) epsilon (double) is the weakening parameter epsilon.
*
* Command line arguments that can be passed in any order:
* -t or --time: Simulation timesteps. Default is 1000. Enter as a simple integer, i.e., 100000 instead of 1e6 or 100,000.
* -s or --size: Simulation size (number of cells). Default is 1000. Enter as a simple integer.
* -r or --rate: Simulation strain (driving) rate. Default value 0.0. Enter as a simple float.
* -d or --disorder: Width of arrest stress distribution. Default value is 0.05. Enter as a simple float.
* -e or --epsilon: Weakening or Strengthening parameter epsilon. Default value is 0.0. Enter as a simple float.
* -o or --output: Type of file to output received as an int. Default is 0 which is both stress and strain, 1 for stress only, 2 for strain only.  
*
* Specific for the Cluster Version
* --taskid: Append the SLURM task ID to the end of output files name
* --jobid: Append the SLURM job ID to the end of output files name
* 
* It write two files as an output:
* - stress_s={-s}_r={-r}_d={-d}_e={-e}.txt is the force on the system and consists of comma-spaced doubles.
* - strain_s={-s}_r={-r}_d={-d}_e={-e}.txt is the strain on the system and consists of comma-spaced doubles.
*/

#include "Slip_Simulation.hpp"

// Mersenne Twister Random Engine Seed 
const static int RAND_SEED = 314159;

// Declares all randomizer objects
static std::mt19937_64 mt_engine{RAND_SEED};
static std::uniform_real_distribution<> uni_rand{0, 1};
static std::weibull_distribution<> wei_rand;

// Declares all global scope variables
static uint32_t time_step = 0;
static double total_stress_at_t = 0;
static double cumulative_total_strain = 0;

// Declares file objects
static std::ofstream stress_file;
static std::ofstream strain_file;

// Variable for storing output file option
static PRINT_TYPE print_option = PRINT_TYPE::BOTH;

/*
* This exists basically only for the sake of origanization.
* The function is inlined into wherever it is called to eliminate call overhead cost.
*/
inline void print_to_file() {

    if (print_option == PRINT_TYPE::STRESS_ONLY) {
        stress_file << std::fixed << std::setprecision(6) << total_stress_at_t << std::endl;
    }
    else if (print_option == PRINT_TYPE::STRAIN_ONLY) {
        strain_file << std::fixed << std::setprecision(6) << cumulative_total_strain << std::endl;
    }
    else {
        stress_file << std::fixed << std::setprecision(6) << total_stress_at_t << std::endl;
        strain_file << std::fixed << std::setprecision(6) << cumulative_total_strain << std::endl;
    }

}

// Main program starts here
int main(int argc, char** argv) {

    constexpr double K = 18.0; // Weibull distribution shape parameter
    const double LAMBDA = 1 / tgamma(1 + 1 / K); // Weibull distribution scale parameter
    constexpr double MODULUS = 0.001; // Elastic shear modulus of the system (J)

    cxxopts::Options options{ "Slip_Simulation", "A simple model simulation of slip-type avalanches." };
    
// Adds commandline options 
#ifdef CLUSTER_BUILD
    options.add_options()
        ("t,time", "Simulation timesteps", cxxopts::value<uint32_t>()->default_value("1000"))
        ("s,size", "Simulation size (num of cells)", cxxopts::value<uint32_t>()->default_value("1000"))
        ("r,rate", "Strain rate", cxxopts::value<double>()->default_value("0.0"))
        ("d,disorder", "Disorder width", cxxopts::value<double>()->default_value("0.05"))
        ("e,epsilon", "Epsilon (Weakening, Strengthing)", cxxopts::value<double>()->default_value("0.0"))
        ("o,output", "Type of file to print out (defaults to BOTH)", cxxopts::value<int>()->default_value("0"))
        ("jobid", "Job (array) ID", cxxopts::value<int>()->default_value("0"))
        ("taskid", "Job array taskid", cxxopts::value<int>()->default_value("0"));
#else
    options.add_options()
        ("t,time", "Simulation timesteps", cxxopts::value<uint32_t>()->default_value("1000"))
        ("s,size", "Simulation size (num of cells)", cxxopts::value<uint32_t>()->default_value("1000"))
        ("r,rate", "Strain rate", cxxopts::value<double>()->default_value("0.0"))
        ("d,disorder", "Disorder width", cxxopts::value<double>()->default_value("0.05"))
        ("e,epsilon", "Epsilon (Weakening, Strengthing)", cxxopts::value<double>()->default_value("0.0"))
        ("o,output", "Type of file to print out (defaults to BOTH)", cxxopts::value<int>()->default_value("0"));
#endif

    auto result = options.parse(argc, argv);

    const uint32_t TIME_MAX = result["time"].as<uint32_t>(); // Max simulation timestep
    const uint32_t AREA = result["size"].as<uint32_t>(); // Total size of the system
    const double RATE = result["rate"].as<double>(); // Driving rate in the case of moving boundary condition
    const double DISORDER = result["disorder"].as<double>(); // Stress redistribution fluctuation rate
    const double EPSILON = result["epsilon"].as<double>(); // Weakening (1 > eps > 0) or Strengthing (-1 < eps < 0) tuning parameter
    const double CONSV = 1 - 1 / sqrt(AREA); // Conservation rate of the system
    
    // Select file type
    if (result["output"].as<int>() == 1) {
        print_option = PRINT_TYPE::STRESS_ONLY;
    }
    else if (result["output"].as<int>() == 2) {
        print_option = PRINT_TYPE::STRAIN_ONLY;
    }
    else {
        print_option = PRINT_TYPE::BOTH;
    }

#ifdef CLUSTER_BUILD
    const int JOBID = result["jobid"].as<int>();
    const int TASKID = result["taskid"].as<int>();
#endif

    // File names
    std::string stress_filename;
    std::string strain_filename;

#ifdef CLUSTER_BUILD

    if (TASKID == 0) {
        stress_filename = "stress_s=" + std::to_string(AREA) + "_r=" + std::to_string(RATE) + "_d=" + std::to_string(DISORDER)
            + "_e=" + std::to_string(EPSILON) + "_jobid=" + std::to_string(JOBID) + ".txt";
        strain_filename = "strain_s=" + std::to_string(AREA) + "_r=" + std::to_string(RATE) + "_d=" + std::to_string(DISORDER)
            + "_e=" + std::to_string(EPSILON) + "_jobid=" + std::to_string(JOBID) + ".txt";
    }
    else {
        stress_filename = "stress_s=" + std::to_string(AREA) + "_r=" + std::to_string(RATE) + "_d=" + std::to_string(DISORDER)
            + "_e=" + std::to_string(EPSILON) + "_jobid=" + std::to_string(JOBID) + "_taskid=" + std::to_string(TASKID) + ".txt";
        strain_filename = "strain_s=" + std::to_string(AREA) + "_r=" + std::to_string(RATE) + "_d=" + std::to_string(DISORDER)
            + "_e=" + std::to_string(EPSILON) + "_jobid=" + std::to_string(JOBID) + "_taskid=" + std::to_string(TASKID) + ".txt";
    }
#else
    stress_filename = "stress_s=" + std::to_string(AREA) + "_r=" + std::to_string(RATE) + "_d=" + std::to_string(DISORDER)
        + "_e=" + std::to_string(EPSILON) + ".txt";
    strain_filename = "strain_s=" + std::to_string(AREA) + "_r=" + std::to_string(RATE) + "_d=" + std::to_string(DISORDER)
        + "_e=" + std::to_string(EPSILON) + ".txt";

#endif

    wei_rand = std::weibull_distribution<>{ K, LAMBDA };

    // Output streams, create file if not existed, replace if existed
    stress_file = std::ofstream{ stress_filename, std::ios::out };
    strain_file = std::ofstream{ strain_filename, std::ios::out };

    /*
    * Technically, these attribute arrays need not be dynamically allocated, but you cannot create normal arrays
    * with non-compile-time-constant size argument (in this case, AREA). This is one of those annoying C++
    * semantics, but it is what it is.
    */
    double* stresses = new double[AREA];
    double* fail_stress = new double[AREA];
    double* arrest_stress = new double[AREA];

    uint32_t i = 0; // Universal index
    bool is_failing = false; // Boolean that is true if the system is in an avalanche, false otherwise

    // Initializing the system
    for (i = 0; i < AREA; i++) {
        arrest_stress[i] = 0.1 * uni_rand(mt_engine) - 0.05;
        fail_stress[i] = 1;
        // This distribution is an approximation of the steady state distribution
        stresses[i] = (1.56585 * pow(i / static_cast<double>(AREA), 0.4) - 0.56585) * (fail_stress[i] - arrest_stress[i]) + arrest_stress[i];
        total_stress_at_t += stresses[i];
    }

    // Main simulation loop
    while (time_step < TIME_MAX) {
        
        double stress_to_fail = 999999; // Minimum stress until the weakest cell fail, initialize to a large number by default
        double redistributed_stress = 0; // Stress to be redistributed
        int failed_count = 0; // Number of cell failed in a timestep

        // If the system is in steady state, "fast forward" to the next avalanche
//-----------------------------------------------------------------------------
        if (!is_failing) {

            // Get minimum stress to fail
            for (i = 0; i < AREA; i++) {
                if (stress_to_fail > fail_stress[i] - stresses[i]) {
                    stress_to_fail = fail_stress[i] - stresses[i];
                }
            }

            // Accumulate stress and strain in non-zero driving rate case
            if (RATE > 0) {
                const uint32_t delta_t = static_cast<uint32_t>(stress_to_fail / (MODULUS * RATE));

                if (delta_t > 0) {
                    const uint32_t t_to_end = (time_step + delta_t) > TIME_MAX - 1 ? TIME_MAX - time_step - 1 : delta_t;
                    
                    // Slowly loads stress and strain until stress reach the critical value
                    for (size_t t = 0; t < t_to_end; t++) {
                        total_stress_at_t += (stress_to_fail / delta_t) * AREA;
                        cumulative_total_strain += RATE;

                        print_to_file();
                    }

                    time_step += t_to_end;
                }

            }
            else { // The case of zero RATE

                cumulative_total_strain += stress_to_fail * MODULUS;

                // Set the system up for an avalanche
                for (i = 0; i < AREA; i++) {
                    stresses[i] += stress_to_fail;
                }
            }
            
        }
//-----------------------------------------------------------------------------

        total_stress_at_t = 0;

        // This is the cell failure and redistribution mechanism
//-----------------------------------------------------------------------------

        // POSSIBLE FUTURE TODO: Add option for uniform random that use DISORDER parameter
 
        // Test for failure and update the appropriate attributes
        for (i = 0; i < AREA; i++) {

            total_stress_at_t += stresses[i];

            if (stresses[i] >= fail_stress[i]) {

                // Stress is lost and redistributed randomly via Weibull distribution
                double lost_stress = (fail_stress[i] - arrest_stress[i]) * wei_rand(mt_engine); 
                cumulative_total_strain += 1 / static_cast<double>(AREA * AREA);
                redistributed_stress += lost_stress;
                stresses[i] -= lost_stress * (1 + CONSV / (AREA - 1));

                if (fail_stress[i] == 1) {
                    fail_stress[i] = 1 - EPSILON * (1 - arrest_stress[i]); // Apply weakening or strengthening in the case EPSILON != 0
                }

                is_failing = true; // Cells failed! Say that the system is in an avalanche
                failed_count++;
            }
        }

        print_to_file();

        time_step++;
        
        // Redistribute stress
        for (i = 0; i < AREA; i++) {
            stresses[i] += redistributed_stress * (CONSV / (AREA - 1)) + (1 / CONSV - 1) * RATE;

            // Say that the avalanche is over when no cell failed during current timestep
            if (failed_count == 0) {
                fail_stress[i] = 1; // Restore weaken or strengthen failure threshold
                is_failing = false;
            }
        }
//-----------------------------------------------------------------------------
        
    }

    // Close files.
    // Not technically neccessary, but since the stream objects are static, might as well
    stress_file.close();
    strain_file.close();

    // De-allocate all dynamic variables.
    // !!!VERY IMPORTANT!!!
    delete[] stresses;
    delete[] fail_stress;
    delete[] arrest_stress;
}

/*
* This is just an archive of the old version's randomizer, included for historical reasons.
* It will NOT be compiled with working code by default, unless the macro _OLD_CODE is defined,
* but frankly, there is no reason to.
*/
#ifdef _OLD_CODE

/*
Copyright (C) 1998 Matthew C. Kuntz, James P. Sethna, Karin A. Dahmen and John Carpenter.
This file is part of the Hysteresis program.
The Hysteresis program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.
See the file COPYING for details.
[yeah I don't have this - Alan]
*/
double* ranmarin(int ijkl, int N) {
    double c, cd, cm, u[97];
    int i97, j97, y;
    double* uni;
    uni = (double*)malloc(N * sizeof(double));
    double output;
    int i, ii, j, jj, k, l, m;
    double s, t;
    int BIGPRIME;
    BIGPRIME = 899999963;
    ijkl = ijkl % BIGPRIME;
    int ij = ijkl / 30082;
    int kl = ijkl - 30082 * ij;

    i = ((ij / 177) % 177) + 2;
    j = (ij % 177) + 2;
    k = ((kl / 169) % 178) + 1;
    l = kl % 169;
    for (ii = 0; ii < 97; ii++) {
        s = 0.0;
        t = 0.5;
        for (jj = 0; jj < 24; jj++) {
            m = (((i * j) % 179) * k) % 179;
            i = j;
            j = k;
            k = m;
            l = (53 * l + 1) % 169;
            if (((l * m) % 64) >= 32) s += t;
            t *= 0.5;
        }
        u[ii] = s;
    }
    c = 362436.0 / 16777216.0;
    cd = 7654321.0 / 16777216.0;
    cm = 16777213.0 / 16777216.0;
    i97 = 96;
    j97 = 32;
    for (y = 0; y < N; y++) {
        uni[y] = u[i97] - u[j97];
        if (uni[y] < 0.0) uni[y] += 1.0;
        u[i97] = uni[y];
        if (--i97 < 0) i97 = 96;
        if (--j97 < 0) j97 = 96;
        c -= cd;
        if (c < 0.0) c += cm;
        uni[y] -= c;
        if (uni[y] < 0.0) uni[y] += 1.0;
    }
    return(uni);
}

/*
This is a makes a Weibull-ly distributed random number.
It uses the uniform rng above.
ijkl is a seed and N is the length of the resultant random array.
k and lambda are the shape parameter and mean respectively.
It outputs uni which is an array of random doubles.
It is based on the code cited below.
*/
double* ranwbl(int ijkl, int N, double k, double lambda) {
    double* uni;
    double* wbl;
    wbl = (double*)malloc(N * sizeof(double));
    int i;
    uni = ranmarin(ijkl, N);
    for (i = 0; i < N; i++) {
        if (uni[i] != 0.0) {
            wbl[i] = lambda * pow(-1.0 * log(uni[i]), 1 / k);
        }
        else {
            wbl[i] = 10.0;
        }
    }
    free(uni);
    return(wbl);
}

#endif