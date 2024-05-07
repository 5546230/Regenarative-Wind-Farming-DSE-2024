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
    dy = r * math.sqrt(3)

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
    height = (row_count * dy) + r  # Extra r for top boundary
    side_length = max(width, height)
    area_of_square = side_length ** 2

    # Create plot
    fig, ax = plt.subplots()
    ax.set_xlim(0, side_length)
    ax.set_ylim(0, side_length)
    ax.set_aspect('equal', 'box')

    # Draw circles
    for (x, y) in positions:
        circle = plt.Circle((x, y), r, edgecolor='blue', facecolor='lightblue', alpha=0.6)
        ax.add_patch(circle)

    #ax.set_title(f"{n} Circles of radius {r} in a staggered pattern in a square of area {area_of_square}")
    ax.set_xlabel("Width [m]")
    ax.set_ylabel("Height [m]")
    plt.show()

    return area_of_square

# Example usage:
n = 10  # Number of circles
r = 2*53.0732*1.05/2   # Radius of each circle
x = int(np.ceil(math.sqrt(n)))   # Base number of circles in the first row
area = draw_and_calculate_variable_hexagonal_square_area(n, r, x)
print(f"The area of the square that contains {n} circles of radius {r} is: {area}")
