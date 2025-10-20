import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import flocksim
rho=0.7
v0=50
eeta=0.12
N=60000
reps=10000
t=0.01
frame_rate=100
def initialisation(N=225, rho=0.509):
    L = np.sqrt(N / rho)
    print("hi")
    print(f"Box side length: {L}")
    coordinates = L * np.random.rand(N, 2) 
    orientation = np.random.rand(N,2)-0.5  
    norm=np.linalg.norm(orientation, axis=1, keepdims=True)
    orientation=orientation/norm
    return coordinates, orientation, L
coordinates,orientation,L=initialisation(N, rho)
coordinates=np.array(coordinates)
orientation=np.array(orientation)
print(type(coordinates), type(orientation))
r_c=L/10
sim=flocksim.Simulation_flock(L, r_c, N)
plot_cords=[]
plot_orien=[]
equilib_orien=None
def integrate(coordinates, orientation,eeta):
    global plot_cords
    global plot_orien  
    global equilib_orien
    for i in range(reps):
        print(i)
        coordinates_n=coordinates+v0*orientation*t
        coordinates=coordinates.ravel().tolist()
        orientation=orientation.ravel().tolist()
        sim.makegrid(coordinates)
        avg_orientation_list=sim.force(coordinates, orientation)
        avg_orientation = np.array(avg_orientation_list).reshape(N, 2)
        avg_angles = np.arctan2(avg_orientation[:, 1], avg_orientation[:, 0])
        random_noise_angles = eeta * np.pi * (2 * np.random.rand(N) - 1)
        new_angles = avg_angles + random_noise_angles
        new_orientation = np.stack([np.cos(new_angles), np.sin(new_angles)], axis=1)
        orientation=new_orientation
        coordinates=np.array(coordinates).reshape(N, 2)
        coordinates=coordinates%L
        if i%frame_rate==0:
             plot_cords.append(coordinates)
             plot_orien.append(orientation)
        coordinates=coordinates_n
    return plot_orien
#integrate(coordinates, orientation)
def calc_order_parameter(plot_orien):
    v_a=[]
    equilib_point=1000  
    for i in plot_orien[equilib_point:]:
        sum_orien=np.sum(i,axis=0)
        avg_vel=v0*sum_orien/N
        v_a.append(np.linalg.norm(avg_vel)/v0)
        avg_va=np.mean(v_a)
    return avg_va
def noise_driven_phase_transition(noise_list,N_list,rho):
    va_list=[]
    for eeta in noise_list:
        for N in N_list:
            orien_list=integrate()

            

def animate(plot_cords, plot_orien, interval=50):
    fig = plt.figure()
    ax = fig.add_subplot()
    def update(frame):
        cords = plot_cords[frame]
        orien = plot_orien[frame]
        ax.clear()
        ax.quiver(cords[:, 0], cords[:, 1], orien[:, 0], orien[:, 1], color='b')
        ax.set_xlim([0, L])
        ax.set_ylim([0, L])
        ax.set_title(f"Time: {frame * frame_rate * t:.2f}s")
        ax.grid(True)
    ani = FuncAnimation(fig, update, frames=len(plot_cords), interval=interval)
    plt.show()
"final_coords, final_orien = integrate(coordinates, orientation)"
animate(plot_cords, plot_orien, 50)

