import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use("TkAgg")
# Compiler optimization levels
import matplotlib.pyplot as plt
import numpy as np

# Optimization levels
# opt_levels = ['-O0', '-O1', '-O2', '-O3']
# x = np.arange(len(opt_levels))

# # Timings for each version
# times = {
#     'v0': [40.4803, 23.5169, 23.8566, 23.5780],
#     'v1': [21.9629, 7.54683, 7.66358, 7.02336],
#     'v2': [15.1869, 5.57040, 5.74219, 4.48152],
#     'v3': [13.9059, 5.56629, 4.37055, 4.33516],
#     'v4': [13.9427, 5.37918, 4.36592, 4.32799],
    
#     'v5': [26.8942, 15.9821, 13.2702, 13.6270],
# }

# # Distinct markers
# markers = {
#     'v0': 'o', 'v1': 's', 'v2': '^',
#     'v3': 'D', 'v4': 'X', 'v5': '*'
# }

# # Slight horizontal offsets for jittering
# offsets = {
#     'v0': 0.0, 'v1': 0.0, 'v2': 0.0,
#     'v3': -0.05, 'v4': 0.05, 'v5': 0.0 
# }

# # Plot setup
# plt.figure(figsize=(10, 6))

# for version, y in times.items():
#     x_shifted = x + offsets[version]
#     plt.plot(x_shifted, y, marker=markers[version], label=version)

# # Axis and legend
# plt.xticks(x, opt_levels)
# plt.xlabel('Compiler Optimization Level')
# plt.ylabel('Time per call (ms)')
# plt.title('compute_forces() Performance')
# plt.legend(title='Version', loc='upper right')
# plt.grid(True)
# plt.tight_layout()

# # Save to PNG
# plt.savefig("compute_forces_performance.png", dpi=300)
# print("Plot saved as compute_forces_performance.png")

# # Show the plot
# plt.show
# import matplotlib.pyplot as plt

# # Data
# threads = [1, 2, 4, 8, 16, 32, 48]
# times = [12718.69, 7965.43, 4084.62, 1739.94, 884.145, 465.149, 319.084]

# # Plot
# plt.figure(figsize=(8, 6))
# plt.plot(threads, times, color='red', marker='o', linestyle='-', linewidth=2, label='Shared Memory Parallel')

# # Add serial performance point
# plt.scatter(1, 4331.080, color='blue', s=100, label='Serial Performance (4331.08s)', zorder=5)

# plt.xlabel('Number of Threads')
# plt.ylabel('Compute Time (s)')
# plt.title('Compute Time vs Number of Threads')
# plt.grid(True)
# plt.xticks(threads)
# plt.legend()
# plt.tight_layout()
# plt.savefig("speed_up.png", dpi=300)
# plt.show()

import matplotlib.pyplot as plt

# Problem sizes
N = [100, 200, 400, 800, 1600, 3200, 6400, 10000]

# CUDA times for different block sizes
cuda_times = {
    16:   [1.976, 3.197, 5.764, 10.923, 20.945, 56.621, 123.855, 230.151],
    32:   [2.216, 3.187, 5.759, 10.954, 21.352, 42.333, 114.452, 185.385],
    64:   [1.923, 3.034, 6.423, 13.595, 27.805, 56.316, 116.586, 192.550],
    128:  [1.994, 3.563, 6.811, 14.712, 31.319, 65.474, 137.903, 218.624],
    256:  [1.937, 4.054, 10.326, 23.058, 48.472, 99.304, 202.911, 321.365]
}

# Serial times
serial_times = [0.439, 1.766, 7.007, 28.305, 111.563, 445.864, 1772.920, 4331.080]

# Plot
plt.figure(figsize=(10, 6))

# Plot CUDA timings
for block_size, times in cuda_times.items():
    plt.plot(N, times, marker='o', linestyle='--', label=f'CUDA Block Size {block_size}')

# Plot serial timings
plt.plot(N, serial_times, marker='o', color='blue', label='Serial')

plt.xlabel('Number of Particles (N)')
plt.ylabel('Compute Time (s)')
plt.title('Runtime vs Number of Particles')
plt.xscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("cuda.png", dpi=300)
plt.show()
