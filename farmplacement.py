#to see how the wind farm locates its multi rotor units

import numpy as np
import matplotlib.pyplot as plt

def plot_wind_farm(num_turbines, area_side_length):
    """
    Plot a map of a wind farm with turbines placed uniformly.

    Args:
    num_turbines (int): Number of turbines in the wind farm.
    area_side_length (float): The length of one side of the square area in km.
    """
    # Calculate the number of turbines per side of the square root of the total number of turbines
    per_side = int(np.ceil(np.sqrt(num_turbines)))
    
    # Generate grid coordinates
    x_coords = np.linspace(0, area_side_length, per_side)
    y_coords = np.linspace(0, area_side_length, per_side)
    
    # Create a meshgrid for turbine positions
    X, Y = np.meshgrid(x_coords, y_coords)
    
    # Flatten the arrays to get the list of coordinates
    x_list = X.ravel()[:num_turbines]
    y_list = Y.ravel()[:num_turbines]
    
    # Plotting
    plt.figure(figsize=(10, 10))
    plt.scatter(x_list, y_list, c='green', marker='o', s=100)  # Turbine marker
    plt.grid(True)
    plt.title('Wind Farm Layout')
    plt.xlabel('Kilometers')
    plt.ylabel('Kilometers')
    plt.axis([0, area_side_length, 0, area_side_length])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

# Example usage
plot_wind_farm(200, 14)  # 200 turbines in a 14km x 14km area
