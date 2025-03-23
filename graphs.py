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
opt_levels = ['-O0', '-O1', '-O2', '-O3']
x = np.arange(len(opt_levels))

# Timings for each version
times = {
    'v0': [40.4803, 23.5169, 23.8566, 23.5780],
    'v1': [21.9629, 7.54683, 7.66358, 7.02336],
    'v2': [15.1869, 5.57040, 5.74219, 4.48152],
    'v3': [13.9059, 5.56629, 4.37055, 4.33516],
    'v4': [13.9427, 5.37918, 4.36592, 4.32799],
    
    'v5': [26.8942, 15.9821, 13.2702, 13.6270],
}

# Distinct markers
markers = {
    'v0': 'o', 'v1': 's', 'v2': '^',
    'v3': 'D', 'v4': 'X', 'v5': '*'
}

# Slight horizontal offsets for jittering
offsets = {
    'v0': 0.0, 'v1': 0.0, 'v2': 0.0,
    'v3': -0.05, 'v4': 0.05, 'v5': 0.0 
}

# Plot setup
plt.figure(figsize=(10, 6))

for version, y in times.items():
    x_shifted = x + offsets[version]
    plt.plot(x_shifted, y, marker=markers[version], label=version)

# Axis and legend
plt.xticks(x, opt_levels)
plt.xlabel('Compiler Optimization Level')
plt.ylabel('Time per call (ms)')
plt.title('compute_forces() Performance')
plt.legend(title='Version', loc='upper right')
plt.grid(True)
plt.tight_layout()

# Save to PNG
plt.savefig("compute_forces_performance.png", dpi=300)
print("Plot saved as compute_forces_performance.png")

# Show the plot
plt.show()
