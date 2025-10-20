import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import flocksim 
class VicsekSimulation:
    def __init__(self,N,L,v0,eeta,N_list=None,eeta_list=None,rho_list=None):
        self.N = N
        self.L= L
        self.v0 = v0
        self.eeta = eeta
        self.rho = self.N / (self.L**2)
        self.r_c = 1  #note that this is the reduced unit of length in this paper
        self.coordinates = self.L * np.random.rand(self.N, 2)
        orientations = np.random.rand(self.N, 2) - 0.5
        norms = np.linalg.norm(orientations, axis=1, keepdims=True)
        self.orientation = orientations / norms
        self.sim_backend = flocksim.Simulation_flock(self.L, self.r_c, self.N)
        self.N_list = N_list
        self.eeta_list = eeta_list
        self.rho_list= rho_list
        self.plot_cords = []
        self.plot_orien = []
        print(f"Initialized simulation with N={self.N}, L={self.L:.2f}")
    def run(self, reps=10000, t=0.01, frame_rate=100):
        for i in range(reps):
            if i % 100 == 0:
                print(f"Step {i}/{reps}")
            coordinates_next = self.coordinates + self.v0 * self.orientation * t
            coords_list = self.coordinates.ravel().tolist()
            orien_list = self.orientation.ravel().tolist()           
            self.sim_backend.makegrid(coords_list)
            avg_orientation_list = self.sim_backend.force(coords_list, orien_list)
            avg_orientation = np.array(avg_orientation_list).reshape(self.N, 2)          
            avg_angles = np.arctan2(avg_orientation[:, 1], avg_orientation[:, 0])
            #noise = self.eeta * np.pi * (2 * np.random.rand(self.N) - 1)
            noise = self.eeta * (np.random.rand(N) - 0.5)
            new_angles = avg_angles + noise
            self.orientation = np.stack([np.cos(new_angles), np.sin(new_angles)], axis=1)
            self.coordinates %= self.L
            if i % frame_rate == 0:
                self.plot_cords.append(self.coordinates.copy())
                self.plot_orien.append(self.orientation.copy())
            self.coordinates = coordinates_next
        print("Simulation finished.")
        return self.plot_orien
    def calculate_order_parameter(self, equilibration_frames=10):
        relevant_frames = self.plot_orien[equilibration_frames:]
        va_values = []
        for snapshot in relevant_frames:
            sum_vector = np.sum(snapshot, axis=0)
            va_instant = np.linalg.norm(sum_vector)/self.N
            va_values.append(va_instant)    
        return np.mean(va_values)
    def noise_driven_phase_transition(self, rho_fixed=4):
        va_list={}
        for N in self.N_list:
            va_list_temp=[]
            L = np.sqrt(N / rho_fixed)
            for eeta in self.eeta_list:
                sim=VicsekSimulation(N=N,L=L,v0=0.03,eeta=eeta)
                sim.run()
                va=sim.calculate_order_parameter(equilibration_frames=10)
                va_list_temp.append(va)
            va_list[N]=va_list_temp#each list inside va_list corresponds to a particular N      
        return va_list
    def density_driven_phase_transition(self,L_fixed=20):
        va_list={}       
        for rho in self.rho_list:
            sim=VicsekSimulation(N=int(round(rho*(L_fixed)**2)),L=L_fixed,v0=0.03,eeta=self.eeta)
            sim.run()
            va=sim.calculate_order_parameter(equilibration_frames=10)
            va_list[rho]=va#each list inside va_list corresponds to a particular rho
        return va_list
    def critical_exponents(self, noise_data, density_data, eta_c_estimates, rho_c_estimate):
        # --- Part 1: Analyze for Beta ---
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)

        for N, va_values in noise_data.items():
            eta_c = eta_c_estimates.get(N)
            if eta_c is None:
                print(f"Warning: No eta_c estimate for N={N}. Skipping.")
                continue

            eeta_list = np.array(self.eeta_list)
            va_values = np.array(va_values)

            # Select data below the critical point
            mask = eeta_list < eta_c
        
            # Calculate the reduced parameter and filter data for log plot
            reduced_eta = (eta_c - eeta_list[mask]) / eta_c
            va_masked = va_values[mask]
        
            log_mask = (reduced_eta > 0) & (va_masked > 0)
            log_x = np.log10(reduced_eta[log_mask])
            log_y = np.log10(va_masked[log_mask])
        
            # Plot the data points
            plt.plot(log_x, log_y, 'o', label=f'N={N}')
        
            # Perform linear regression to find the slope (beta)
            if len(log_x) > 1:
                slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)
                print(f"For N={N}, Beta (slope) = {slope:.3f}")

        plt.xlabel('log [ ($\eta_c - \eta) / \eta_c$ ]')
        plt.ylabel('log [ $v_a$ ]')
        plt.title('Critical Exponent $\\beta$')
        plt.legend()
        plt.grid(True)

    # --- Part 2: Analyze for Delta ---
        plt.subplot(1, 2, 2)
    
        rho_list = np.array(density_data['density'])
        va_values = np.array(density_data['va'])

        # Select data above the critical point
        mask = rho_list > rho_c_estimate
    
        # Calculate the reduced parameter and filter data
        reduced_rho = (rho_list[mask] - rho_c_estimate) / rho_c_estimate
        va_masked = va_values[mask]
    
        log_mask = (reduced_rho > 0) & (va_masked > 0)
        log_x = np.log10(reduced_rho[log_mask])
        log_y = np.log10(va_masked[log_mask])

    # Perform linear regression to find the slope (delta)
        if len(log_x) > 1:
            slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)
            print(f"For density transition, Delta (slope) = {slope:.3f}")
            # Plot the fitted line
            fit_line = slope * log_x + intercept
            plt.plot(log_x, fit_line, 'r-', label=f'Fit (slope={slope:.2f})')
        
        plt.plot(log_x, log_y, 's')
        plt.xlabel('log [ ($\\rho - \\rho_c) / \\rho_c$ ]')
        plt.ylabel('log [ $v_a$ ]')
        plt.title('Critical Exponent $\\delta$')
        plt.legend()
        plt.grid(True)
    
        plt.tight_layout()
        plt.show()

    def animate(self, t=0.01, frame_rate=100, interval=50):
        fig, ax = plt.subplots()
        def update(frame):
            cords = self.plot_cords[frame]
            orien = self.plot_orien[frame]
            ax.clear()
            ax.quiver(cords[:, 0], cords[:, 1], orien[:, 0], orien[:, 1], color='b', headwidth=2, headlength=3)
            ax.set_xlim([0, self.L])
            ax.set_ylim([0, self.L])
            time = frame * frame_rate * t
            ax.set_title(f"Time: {time:.2f}s | $v_a$: {self.calculate_order_parameter():.3f}")
            ax.grid(True)       
        ani = FuncAnimation(fig, update, frames=len(self.plot_cords), interval=interval)
        plt.show()
if __name__ == "__main__":
    N = 5000
    L = 25
    v0 = 0.03
    eeta = 0.1
    N_list=[40,100,400,4000,10000]
    eeta_list=np.linspace(0,1,50)
    rho_list=np.linspace(0,10,30)
    simulation = VicsekSimulation(N=N,L=L,v0=v0, eeta=eeta, N_list=N_list, eeta_list=eeta_list, rho_list=rho_list)
    #simulation.run(reps=5000, frame_rate=50)
    va_noise = simulation.noise_driven_phase_transition()
    va_density = simulation.density_driven_phase_transition(L_fixed=20)
    #plotting va vs eeta for different N
    plt.figure(figsize=(8,5))
    for N, va_list in va_noise.items():
        plt.plot(simulation.eeta_list, va_list, marker='o', linestyle='-', label=f"N={N}")
    plt.xlabel("η (noise)")
    plt.ylabel("<v_a>")
    plt.title("Order parameter vs noise")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plotting va vs rho
    rho_items = sorted(va_density.items(), key=lambda x: x[0])
    rhos = [item[0] for item in rho_items]
    vas = [item[1] for item in rho_items]
    plt.figure(figsize=(8,5))
    plt.plot(rhos, vas, marker='o', linestyle='-')
    plt.xlabel("Density ρ")
    plt.ylabel("<v_a>")
    plt.title("Order parameter vs density")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    order_param = simulation.calculate_order_parameter(equilibration_frames=10)
    print(f"\nFinal time-averaged order parameter <v_a>: {order_param:.4f}")
    simulation.animate(frame_rate=200)