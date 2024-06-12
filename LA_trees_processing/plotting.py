import numpy as np
import matplotlib.pyplot as plt

# Function to read data from a file
def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1]

# File names
files = ['c_lorentzFit.txt', 'c_lorentzFit380947_L1B.txt', 'c_lorentzFit380947_L1F.txt']
labels = ['Lorentz Fit', 'Lorentz Fit 380947 L1B', 'Lorentz Fit 380947 L1F']

# Initialize plot
plt.figure(figsize=(10, 6))

# Read and plot data from each file
for file, label in zip(files, labels):
    x, y = read_data(file)
    plt.plot(x, y, label=label)

# Customize plot
plt.xlabel('X values')
plt.ylabel('Y values')
plt.title('Plot of X vs Y from Three Files')
plt.legend()
plt.grid(True)

# Show plot
plt.show()
