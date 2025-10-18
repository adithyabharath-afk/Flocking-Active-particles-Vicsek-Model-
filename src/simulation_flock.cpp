#include "simulation_flock.h" // Include the header to ensure consistency
#include <cmath>
#include <vector>
#include <iostream> // Good for debugging if needed

// --- Important ---
// The class definition is GONE from this file. We are only IMPLEMENTING
// the functions that were DECLARED in simulation.h.

// Use the scope resolution operator '::' to tell the compiler these
// functions belong to the 'simulation' class.
// Constructor implementation. Note: num_particles is now an 'int' to match the header.
simulation_flock::simulation_flock(double Length, double cutoff, int num_particles) {
    L = Length;
    r_c = cutoff;
    N = num_particles;
    // It's good practice to initialize all member variables
    nx = 0;
    ny = 0;
}

// makegrid implementation. Note: 'coordinates' is now 'const' to match the header.
void simulation_flock::makegrid(const std::vector<double>& coordinates) {
    nx = static_cast<int>(floor(L / r_c));
    ny = nx;
    grid.assign(nx, std::vector<std::vector<int>>(ny));

    for (int i = 0; i < N; ++i) {
        double x = coordinates[2 * i];
        double y = coordinates[2 * i + 1];
        int cellx = static_cast<int>(floor(x / r_c));
        int celly = static_cast<int>(floor(y / r_c));
        if (cellx >= 0 && cellx < nx && celly >= 0 && celly < ny) {
            grid[cellx][celly].push_back(i);
        }
    }
}

// force2dhp implementation. Note: 'coordinates' is now 'const' to match the header.
//sphagetti code
std::vector<double> simulation_flock::force(const std::vector<double>& coordinates,const std::vector<double>& orientations) {
    double orientationx=0.0,orientationy=0.0, r_ij_x, r_ij_y;
    std::vector<double> particles_in_range(N,0.0);
    double r_c2=r_c*r_c;
    std::vector<double> avg_orientation(2*N,0.0);
    for (int i=0;i<nx;i++){
        for (int j=0;j<ny;j++){
            for (int dx=-1;dx<=1;dx++){
                for (int dy=-1;dy<=1;dy++){
                    int ln_x = (i + dx + nx) % nx;//ln:local neighbour
                    int ln_y = (j + dy + ny) % ny;    
                    for(int k:grid[i][j]){
                        for(int l:grid[ln_x][ln_y]){
                            r_ij_x = coordinates[2 * k] - coordinates[2 * l];
                            r_ij_y = coordinates[2 * k + 1] - coordinates[2 * l + 1];
                            r_ij_x -= L * round(r_ij_x / L);
                            r_ij_y -= L * round(r_ij_y / L);
                            double r_sq = r_ij_x * r_ij_x + r_ij_y * r_ij_y;
                            if (r_sq > r_c2) continue;
                            if (k==l) continue;
                            particles_in_range[k]+=1;
                            orientationx+=orientations[2*l];
                            orientationy+=orientations[2*l+1];
                        }
                        avg_orientation[2*k]+=orientationx;
                        avg_orientation[2*k+1]+=orientationy;                      
                        orientationx=0.0;
                        orientationy=0.0;
                    }    
                }
            }
            for (int k:grid[i][j]){
                if(particles_in_range[k]>0){
                    avg_orientation[2*k]=avg_orientation[2*k]/particles_in_range[k];
                    avg_orientation[2*k+1]=avg_orientation[2*k+1]/particles_in_range[k];
                }
                //particles_in_range[k]=0.0;
            }
        }
    }
    return avg_orientation;
};

