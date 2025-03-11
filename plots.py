import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
matplotlib.use("TkAgg")  # Ensure plots open in a separate window

# Load kinetic energy data and plot
def plot_kinetic_energy(filename="kinetic_energy.txt"):
    data = np.loadtxt(filename)
    time = data[:, 0]
    kinetic_energy = data[:, 1]

    plt.figure(figsize=(8, 5))
    plt.plot(time, kinetic_energy, marker="o", linestyle="-", color="b", markersize=3)
    plt.xlabel("Time (units)")
    plt.ylabel("Kinetic Energy")
    plt.title("Kinetic Energy vs. Time")
    plt.grid(True)
    plt.show(block=True)  # Ensure plot opens in a new window

# Load particle data and create animation
def animate_xy_plane(filename="particles.txt"):
    data = np.loadtxt(filename)

    # Get unique time steps
    unique_times = np.unique(data[:, 0])

    # Create figure for animation
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(np.min(data[:, 2]), np.max(data[:, 2]))  # X limits
    ax.set_ylim(np.min(data[:, 3]), np.max(data[:, 3]))  # Y limits
    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_title("Particle Positions in XY Plane (Animation)")

    # Initialize scatter plot with different colors
    scatter = ax.scatter([], [], c=[])

    # Function to update frame
    def update(frame):
        time_step = unique_times[frame]
        mask = data[:, 0] == time_step
        x = data[mask, 2]  # X positions
        y = data[mask, 3]  # Y positions
        types = data[mask, 1]  # Particle types (assumed stored in column 1)

        # Assign colors based on particle type
        colors = np.where(types == 0, "blue", "red")

        scatter.set_offsets(np.c_[x, y])
        scatter.set_color(colors)  # Update colors

        ax.set_title(f"Particle Positions at t = {time_step:.2f}")
        return scatter,

    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(unique_times), interval=100, blit=False)

    plt.show(block=True)  # Ensure animation opens in a new window

# Run plots
if __name__ == "__main__":
    plot_kinetic_energy()
    animate_xy_plane()
