import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import flocksim
print("hi")
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
def integrate(coordinates, orientation):
    global plot_cords
    global plot_orien  
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
integrate(coordinates, orientation)
'''def animate(plot_cords,interval=5):
    fig=plt.figure()#creating a blank canvas
    ax=fig.add_subplot()#creates a specific plotting area in the canvas
    def update(frame):
        cords=plot_cords[frame]       
        ax.clear()#erases everything on the ax plotting area,the points,the tilte , the dimentions everything       
        ax.scatter(cords[:,0],cords[:,1],s=1)
        ax.set_xlim([0,L])
        ax.set_ylim([0,L])
        ax.set_title(f"time{frame*frame_rate:.2f}")
        ax.grid(True)
    ani=FuncAnimation(fig,update,frames=len(plot_cords),interval=50)
    plt.show()
animate(plot_cords,50)'''
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
