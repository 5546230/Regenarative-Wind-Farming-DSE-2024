import matplotlib.pyplot as plt
import math
import numpy as np

def calculate_horizontal_overlap_distance(r, overlap_area):
    """Calculate the horizontal distance for two circles to achieve a specific overlap area."""
    def intersection_area(d):
        return 2 * r**2 * math.acos(d / (2 * r)) - 0.5 * d * math.sqrt(4 * r**2 - d**2)

    target_area = overlap_area * math.pi * r**2
    low, high = 0, 2 * r
    while high - low > 1e-5:
        mid = (low + high) / 2
        if intersection_area(mid) < target_area:
            high = mid
        else:
            low = mid
    return mid

def draw_and_calculate_hexagonal_area_no_vertical_overlap(n, r, x):
    """Draws 'n' circles with horizontal overlap and calculates the square area."""
    # Initialize variables to track row details and total circles placed
    num_circles_placed = 0
    row_count = 0
    positions = []

    # Calculate the horizontal overlap distance for a 26% area overlap
    horizontal_overlap_distance = calculate_horizontal_overlap_distance(r, 0.26)

    # Calculate horizontal distance based on horizontal overlap and vertical distance as twice the radius
    dx = horizontal_overlap_distance
    dy = 1.9 * r  # no vertical overlap, just touch each other vertically

    # Generate positions of circles row by row
    while num_circles_placed < n:
        row_length = x if row_count % 2 == 0 else x + 1
        if num_circles_placed + row_length > n:
            row_length = n - num_circles_placed  # Adjust last row if it exceeds total circles

        # Calculate offset for centering staggered rows
        offset = 0 if row_count % 2 == 0 else dx / 2
        
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
        circle = plt.Circle((x, y), r, facecolor='lightblue', alpha=0.6)
        ax.add_patch(circle)

    ax.set_xlabel("Width [m]")
    ax.set_ylabel("Height [m]")
    plt.show()

    return area_of_square

# Example usage:
n = 50  # Number of circles
r = 23.13406979*1.05   # Radius of each circle
x = int(np.ceil(np.sqrt(n)))   # Base number of circles in the first row
area = draw_and_calculate_hexagonal_area_no_vertical_overlap(n, r, x)
print(f"The area of the square that contains {n} circles of radius {r} is: {area}")

