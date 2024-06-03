import math
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import csv
import os
from typing import List, Dict, Any

def custom_map():
    colors = [
        (0, 'blue'),
        (0.5, 'darkgray'),
        (1, 'red')
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


class CsvOutput:
    def __init__(self, fname: str, output_dir: str = '../Outputs/'):
        self.fname = fname
        self.filepath = os.path.abspath(output_dir+fname)

    def write(self, data: Dict[str, List[Any]]):
        with open(self.filepath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            header = ['Rod'] + list(data.keys())
            writer.writerow(header)
            num_rods = len(next(iter(data.values())))

            for i in range(num_rods):
                row = [i + 1] + [data[param][i] for param in data]
                writer.writerow(row)

    def read(self):
        data = {}
        with open(self.filepath, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)

            header = next(reader)
            param_names = header[1:]

            for param in param_names:
                data[param] = []

            for row in reader:
                for i, param in enumerate(param_names):
                    data[param].append(float(row[i + 1]))
        return data

    def check_current(self)->bool:
        return os.path.exists(self.filepath) and os.path.getsize(self.filepath) > 0



if __name__ == "__main__":
    fem_output = {
        'stress': [100, 200, 150, 800],
        'strain': [0.01, 0.02, 0.015, .1],
        'displacement': [0.1, 0.2, 0.15, .1]
    }

    csv_output = CsvOutput(fname=r'fem_results_1.csv', output_dir = '../Outputs/')
    #csv_output.write(fem_output)

    if csv_output.check_current():
        read_data = csv_output.read()
        print(read_data)