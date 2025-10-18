#ifndef SIMULATION_FLOCK_H
#define SIMULATION_FLOCK_H

#include <vector>
#include <utility> // Required for std::pair

// Declaration of the simulation class. This is the "menu".
class simulation_flock {
public:
    // Member variables
    double L, r_c;
    int nx, ny, N;
    std::vector<std::vector<std::vector<int>>> grid;

    // Constructor declaration
    simulation_flock(double Length, double cutoff, int num_particles);

    // Method declarations
    void makegrid(const std::vector<double>& coordinates);
    std::vector<double> force(const std::vector<double>& coordinates,const std::vector<double>& orientations);
};

#endif // SIMULATION_FLOCK_H