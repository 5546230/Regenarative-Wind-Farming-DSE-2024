import matplotlib.pyplot as plt
import math
import numpy as np

def draw_and_calculate_variable_hexagonal_square_area(n, r, x):
    """
    Draws 'n' circles of radius 'r' with a row alternating configuration between X and X+1 circles
    and calculates the area of the minimum square that fully contains them.
    
    :param n: Number of circles
    :param r: Radius of each circle
    :param x: Base number of circles in the first row
    :return: Area of the square that fully contains all the circles
    """
    # Initialize variables to track row details and total circles placed
    num_circles_placed = 0
    row_count = 0
    positions = []

    # Calculate horizontal and vertical distance between centers
    dx = 2 * r
    dy = r *math.sqrt(3) #bc height of equilateral triangle is sqrt(3)/2*2r
    # Generate positions of circles row by row
    while num_circles_placed < n:
        row_length = x if row_count % 2 == 0 else x - 1
        if num_circles_placed + row_length > n:
            row_length = n - num_circles_placed  # Adjust last row if it exceeds total circles

        # Calculate offset for centering and staggered rows
        offset = 0 if row_count % 2 == 0 else r
        
        for i in range(row_length):
            positions.append((dx * i + offset + r, dy * row_count + r))
        
        num_circles_placed += row_length
        row_count += 1

    # Calculate the required dimensions
    width = max(pos[0] for pos in positions) + r  # Extra r for boundary
    height = max(pos[1] for pos in positions) + r  # Extra r for top boundary
    side_length = max(width, height)
    area_of_square = width * height #side_length ** 2
    print(width, height)
    # Create plot
    fig, ax = plt.subplots()
    ax.set_xlim(0, side_length)
    ax.set_ylim(0, side_length)
    ax.set_aspect('equal', 'box')

    # Draw circles
    for (x, y) in positions:
        circle = plt.Circle((x, y), r,  facecolor='lightblue', alpha=0.6)
        rotor = plt.Circle((x, y), r/1.05, edgecolor='red', facecolor='lightblue', alpha=0.6)
        ax.add_patch(circle)
        ax.add_patch(rotor)

    #ax.set_title(f"{n} Circles of radius {r} in a staggered pattern in a square of area {area_of_square}")
    ax.set_xlabel("Width [m]")
    ax.set_ylabel("Height [m]")
    plt.show()

    return area_of_square, positions

# Example usage:
n = 32  # Number of circles
r1 = 170 #radius of single rotor
pg = 0 #power gain due to multi-rotors
r = r1/((n*(1+pg))**0.5)*1.05  # Radius of each circle
x = int(np.ceil(np.sqrt(n)))   # Base number of circles in the first row
area, positions = draw_and_calculate_variable_hexagonal_square_area(n, r, x)
radius1 = r/1.05
print(f"The area of the square that contains {n} circles of radius {radius1} is: {area}")
positions = np.array(positions)
print(positions[8,0])
positions[:,0] = positions[:,0]-positions[8,0]
plt.plot(positions[:,0], positions[:,1], marker = "o")
plt.show()
d = 20
def calculate_I_z_array(x_pos, d):
    
    m_per_turbine = 1280085.241/33+ (2235959.595 / 33)
    I_z = np.sum(m_per_turbine * (x_pos**2))+ m_per_turbine*33*(d**2)

    return I_z
print("Mass moment of inertia")
print(calculate_I_z_array(positions[:,0], d))