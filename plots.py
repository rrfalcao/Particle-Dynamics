
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

matplotlib.use("TkAgg")  # Ensure plots open in a separate window
limit =20
# Load kinetic energy data and plot
def plot_kinetic_energy(filename="kinetic_energy.txt", save_path="kinetic_energy.png"):
    data = np.loadtxt(filename)
    time = data[:, 0]
    kinetic_energy = data[:, 1]

    plt.figure(figsize=(8, 5))
    plt.plot(time, kinetic_energy, marker="o", linestyle="-", color="b", markersize=3)
    plt.xlabel("Time (units)")
    plt.ylabel("Kinetic Energy")
    plt.title("Kinetic Energy vs. Time")
    plt.grid(True)
    plt.tight_layout()
    
    # Save and show
    plt.savefig(save_path)
    plt.show(block=True)

# Load particle data and create animation
def animate_xy_plane(filename="particles.txt"):
    data = np.loadtxt(filename)

    # Get unique time steps
    unique_times = np.unique(data[:, 0])

    # Create figure for animation
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0,limit)  # X limits
    ax.set_ylim(0,limit)  # Y limits
    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_title("Particle Positions in XY Plane (Animation)")

    # Initialize scatter plot
    scatter = ax.scatter([], [], c=[], marker='o')

    # Function to update frame
    def update(frame):
        time_step = unique_times[frame]
        mask = data[:, 0] == time_step
        x = data[mask, 3]  # X positions (column 3)
        y = data[mask, 4]  # Y positions (column 4)
        types = data[mask, 2]  # **Correct column for particle type (column 2)**

        # Assign colors based on particle type
        colors = np.where(types == 0, "blue", "red")

        scatter.set_offsets(np.c_[x, y])
        scatter.set_color(colors)  # Update colors

        ax.set_title(f"Particle Positions at t = {time_step:.2f}")
        return scatter,

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(unique_times), interval=100, blit=False)

    plt.show(block=True)  # Ensure animation opens in a new window


# 3D Particle Trajectory Animation
def animate_3d_trajectory(filename="particles.txt"):
    data = np.loadtxt(filename)

    # Get unique time steps
    unique_times = np.unique(data[:, 0])

    # Create 3D figure
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_zlabel("Z Position")
    ax.set_title("Particle Positions in 3D (Animation)")

    # Set axis limits based on data
    ax.set_xlim(0,limit)  # X limits
    ax.set_ylim(0,limit)  # Y limits
    ax.set_zlim(0,limit)  # Z limits

    # Initialize scatter plot for particles
    scatter = ax.scatter([], [], [], c=[], marker='o')

    # Function to update frame
    def update(frame):
        time_step = unique_times[frame]
        mask = data[:, 0] == time_step
        x = data[mask, 3]  # X positions
        y = data[mask, 4]  # Y positions
        z = data[mask, 5]  # Z positions
        types = data[mask, 2]  # Particle types (column 2)

        # Assign colors based on particle type
        colors = np.where(types == 0, "blue", "red")

        scatter._offsets3d = (x, y, z)  # Update particle positions
        scatter.set_color(colors)  # Update colors
        ax.set_title(f"Particle Positions at t = {time_step:.2f}")

        return scatter,

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(unique_times), interval=100, blit=False)

    plt.show(block=True)  # Ensure animation opens in a new window
# Static 2D XY trajectory plot
def plot_xy_trajectories(filename="particles.txt", save_path="xy_trajectories.png"):
    data = np.loadtxt(filename)
    unique_ids = np.unique(data[:, 1])  # Particle IDs (assumed column 1)

    plt.figure(figsize=(8, 6))

    for i, pid in enumerate(unique_ids):
        mask = data[:, 1] == pid
        x = data[mask, 3]  # X positions
        y = data[mask, 4]  # Y positions

        # Choose color
        if i == 0:
            edgecolor = 'blue'
            label = f"Particle {int(pid)}"
        elif i == 1:
            edgecolor = 'red'
            label = f"Particle {int(pid)}"
        else:
            continue  # Skip others for clarity

        # Plot hollow circles
        plt.scatter(x, y, s=40, facecolors='none', edgecolors=edgecolor, label=label, linewidths=1.5)

    plt.xlabel("X Position")
    plt.ylabel("Y Position")
    plt.title("2D Particle Trajectories - Test 3")
    plt.grid(True)
    plt.legend(loc='upper right', fontsize='medium')
    plt.xlim(0, limit)
    plt.ylim(0, limit)
    plt.tight_layout()

    plt.savefig(save_path)
    plt.show(block=True)
# Run plots
if __name__ == "__main__":
    plot_kinetic_energy()
    plot_xy_trajectories()
    animate_xy_plane()
    # animate_3d_trajectory()
