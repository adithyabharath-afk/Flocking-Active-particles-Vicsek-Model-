#include "simulation_flock.h" // Include the header to ensure consistency
#include <cmath>
#include <vector>
#include <iostream> 
#include <omp.h> 

simulation_flock::simulation_flock(double Length, double cutoff, int num_particles) {
    L = Length;
    r_c = cutoff;
    N = num_particles;
    nx = 0;
    ny = 0;
}
void simulation_flock::makegrid(const std::vector<double>& coordinates) {
    nx = static_cast<int>(floor(L / r_c));
    ny = nx;
    grid.assign(nx, std::vector<std::vector<int>>(ny));
    for (int i = 0; i < N; ++i) {
        double x = coordinates[2 * i];
        double y = coordinates[2 * i + 1];
        int cellx = static_cast<int>(x / r_c);
        int celly = static_cast<int>(y / r_c);
        if (cellx >= nx) cellx = nx - 1;
        if (celly >= ny) celly = ny - 1;
        if (cellx < 0) cellx = 0; 
        if (celly < 0) celly = 0; 
        grid[cellx][celly].push_back(i);      
    }
}
std::vector<double> simulation_flock::force(const std::vector<double>& coordinates,const std::vector<double>& orientations) {
    double orientationx=0.0,orientationy=0.0, r_ij_x, r_ij_y;
    std::vector<double> particles_in_range(N,0.0);
    double r_c2=r_c*r_c;
    std::vector<double> avg_orientation(2*N,0.0);
    #pragma omp parallel for collapse(2)
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
/*std::vector<double> simulation_flock::force(const std::vector<double>& coordinates, const std::vector<double>& orientations) {
    // Vector to store the final average orientation for each particle.
    std::vector<double> avg_orientation(2 * N, 0.0);
    const double r_c2 = r_c * r_c; // Pre-calculate the squared cutoff distance.

    // --- The MAIN LOOP: Iterate "particle-by-particle" ---
    for (int k = 0; k < N; ++k) {

        // --- 1. SETUP: Initialize sums and counters fresh for EACH particle 'k' ---
        double sum_vx = 0.0;
        double sum_vy = 0.0;
        int neighbor_count = 0;

        // Store particle k's coordinates in local variables for faster access.
        const double kx = coordinates[2 * k];
        const double ky = coordinates[2 * k + 1];

        // Find which cell particle 'k' is in.
        const int cellx = static_cast<int>(floor(kx / r_c));
        const int celly = static_cast<int>(floor(ky / r_c));

        // --- 2. WORK: Search the 3x3 neighborhood to find all true neighbors ---
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                
                // Get the coordinates of the neighbor cell, handling periodic boundaries.
                int neighbor_cell_x = (cellx + dx + nx) % nx;
                int neighbor_cell_y = (celly + dy + ny) % ny;

                // Iterate through all candidate particles 'l' in this neighbor cell.
                for (int l : grid[neighbor_cell_x][neighbor_cell_y]) {
                    
                    // Calculate distance with periodic boundary conditions (minimum image convention).
                    double r_ij_x = kx - coordinates[2 * l];
                    double r_ij_y = ky - coordinates[2 * l + 1];
                    r_ij_x -= L * round(r_ij_x / L);
                    r_ij_y -= L * round(r_ij_y / L);
                    
                    double r_sq = r_ij_x * r_ij_x + r_ij_y * r_ij_y;

                    // If 'l' is a true neighbor (within the cutoff radius)...
                    if (r_sq < r_c2) {
                        // ...add its orientation to our sum and increment the count.
                        sum_vx += orientations[2 * l];
                        sum_vy += orientations[2 * l + 1];
                        neighbor_count++;
                    }
                }
            }
        }

        // --- 3. FINISH: After checking all neighbors, perform ONE division ---
        if (neighbor_count > 0) {
            avg_orientation[2 * k] = sum_vx / neighbor_count;
            avg_orientation[2 * k + 1] = sum_vy / neighbor_count;
        }
        // If a particle has no neighbors, its avg_orientation remains (0,0) by default.
    }

    return avg_orientation;
}*/

