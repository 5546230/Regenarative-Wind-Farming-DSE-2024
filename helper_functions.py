import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import math
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

def custom_map():
    colors = [
        (0, 'blue'),  # Blue for negative values
        (0.5, 'darkgray'),  # Gray for zero
        (1, 'red')  # Red for positive values
    ]
    cmap = LinearSegmentedColormap.from_list('CustomMap', colors)
    return cmap


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
            positions.append([dx * i + offset + r, dy * row_count + r])
        num_circles_placed += row_length
        row_count += 1
    width = max(pos[0] for pos in positions) + r
    height = max(pos[1] for pos in positions) + r
    area = width * height
    return positions, width, height, area