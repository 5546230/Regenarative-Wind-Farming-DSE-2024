import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
def custom_map():
    colors = [
        (0, 'blue'),  # Blue for negative values
        (0.5, 'darkgray'),  # Gray for zero
        (1, 'red')  # Red for positive values
    ]
    cmap = LinearSegmentedColormap.from_list('CustomMap', colors)
    return cmap