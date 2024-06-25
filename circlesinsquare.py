import matplotlib.pyplot as plt
import math
import numpy as np

def calculate_grid_positions(n, r, x):
    """Calculate positions for circles arranged in a grid."""
    positions = []
    row_count = 0
    num_circles_placed = 0
    while num_circles_placed < n:
        row_length = x
        if num_circles_placed + row_length > n:
            row_length = n - num_circles_placed
        for i in range(row_length):
            positions.append((2 * r * i + r, 2 * r * row_count + r))
        num_circles_placed += row_length
        row_count += 1
    width = max(pos[0] for pos in positions) + r
    height = max(pos[1] for pos in positions) + r
    area = width * height
    return positions, width, height, area

def calculate_hexagonal_positions(n, r, x):
    """Calculate positions for circles arranged hexagonally."""
    positions = []
    row_count = 0
    num_circles_placed = 0
    dx = 2 * r
    dy = r * math.sqrt(3)
    while num_circles_placed < n:
        row_length = x if row_count % 2 == 0 else x -1
        if num_circles_placed + row_length > n:
            row_length = n - num_circles_placed
        offset = 0 if row_count % 2 == 0 else r
        for i in range(row_length):
            positions.append((dx * i + offset + r, dy * row_count + r))
        num_circles_placed += row_length
        row_count += 1
    width = max(pos[0] for pos in positions) + r
    height = max(pos[1] for pos in positions) + r
    area = width * height
    return positions, width, height, area

def plot_circles(positions, width, height, title):
    """Plot circles with given positions."""
    fig, ax = plt.subplots()
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect('equal', 'box')
    for (x, y) in positions:
        circle = plt.Circle((x, y), r, facecolor='lightblue', alpha=0.6)
        circle2 = plt.Circle((x, y), r/1.05, facecolor='red', alpha=0.6)
        ax.add_patch(circle)
        ax.add_patch(circle2)
    ax.set_title(title)
    ax.set_xlabel("Width [m]")
    ax.set_ylabel("Height [m]")
    plt.show()

# Example usage:
n = 33  # Number of circles
r1 = 170  # Radius of single rotor
pg = 0   # Power gain due to multi-rotors
r = 20.1 #r1 / ((n * (1 + pg)) ** 0.5) * 1.05  # Effective radius of each circle
x = int(np.ceil(np.sqrt(n)))   # Base number of circles in the first row

# Calculate positions and areas
grid_positions, grid_width, grid_height, grid_area = calculate_grid_positions(n, r, x)
hex_positions, hex_width, hex_height, hex_area = calculate_hexagonal_positions(n, r, x)

# Calculate area efficiency
radius1 = r / 1.05
grid_area_efficiency = np.pi * radius1 ** 2 * n / grid_area
hex_area_efficiency = np.pi * radius1 ** 2 * n / hex_area

print(f"Grid layout: Area = {grid_area}, Efficiency = {grid_area_efficiency}")
print(f"Hexagonal layout: Area = {hex_area}, Efficiency = {hex_area_efficiency}")
print("Radius = ", radius1)

# Plotting both layouts
plot_circles(grid_positions, grid_width, grid_height, "Grid Layout")
plot_circles(hex_positions, hex_width, hex_height, "Hexagonal Layout")
print(hex_width, hex_height)

#import csv

# The provided array
data = hex_positions
#print(data)
'''
# The name of the CSV file
csv_file = 'hexpositions_data.csv'

# Writing the data to a CSV file
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y'])  # Writing the header
    writer.writerows(data)  # Writing the data

#csv_file
'''

#===================ESTIMATE I_zz===========================

I_zz_truss = 3.16e+10

m_RNA_single = 16650

coordinates = np.array([[ 30.38 ,       30.38      ],
 [ 30.38      ,  91.14      ],
 [ 30.38      , 151.9       ],
 [ 30.38      , 212.66      ],
 [ 30.38      , 273.42      ],
 [ 30.38      , 334.18      ],
 [ 82.99970353,  60.76      ],
 [ 82.99970353, 121.52      ],
 [ 82.99970353, 182.28      ],
 [ 82.99970353, 243.04      ],
 [ 82.99970353, 303.8       ],
 [135.61940707,  30.38      ],
 [135.61940707,  91.14      ],
 [135.61940707, 151.9       ],
 [135.61940707, 212.66      ],
 [135.61940707, 273.42      ],
 [135.61940707, 334.18 ]])

I_zz_RNA = np.sum(m_RNA_single*(coordinates[:, 0]**2))
m_RNA_total = 34*m_RNA_single
print(I_zz_RNA)
I_zz_AFC = 205e3/12*(277.24**2+50**2)
m_AFC = 205e3
print(I_zz_AFC)

m_tower = 3.033e6 
I_zz_tower = 2.893e+07
I_zz_tot = I_zz_truss+I_zz_AFC+I_zz_RNA#+I_zz_tower
m_tot = m_RNA_total+ m_AFC+m_tower+4315.091e3 +1.765e6

print(f'Total mass moment of inertia: {I_zz_tot}')
print(f'Total mass: {m_tot}')


