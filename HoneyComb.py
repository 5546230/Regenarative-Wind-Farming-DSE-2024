import numpy as np
import pickle
import matplotlib.pyplot as plt

# Load the interpolated function
with open('interpolated_function.pkl', 'rb') as f:
    CTfunction = pickle.load(f)

class Rotor:
    def __init__(self, r, position, grid, q, pitch_angle=0.0):
        self.A = np.pi * r ** 2
        self.position = position
        self.q = q
        self.grid = grid
        self._pitch_angle = pitch_angle
        self.thrust = 0.0
        self.moments = (0.0, 0.0)
        self.update_properties()

    @property
    def pitch_angle(self):
        return self._pitch_angle
    
    @pitch_angle.setter
    def pitch_angle(self, value):
        if self._pitch_angle != value:
            self._pitch_angle = value
            self.update_properties()

    def update_properties(self):
        # Update dependent properties based on the pitch angle
        self.thrust = self.calculate_thrust()
        self.moments = self.calculate_moments()

    def calculate_thrust(self):
        C_T = CTfunction(self._pitch_angle)
        return self.q * C_T * self.A * np.cos(self.grid.yaw_angle) ** 2

    def calculate_moments(self):
        x, z = self.position
        return (self.thrust * z, self.thrust * x)


class HoneycombGrid:
    def __init__(self, config):
        self.n = config['n']
        self.r = config['r']
        self.x = config['x']
        self.q = config['q']
        self.W = config['W']
        self.J = config['J']
        self.R = config['R']
        self.mu = config['mu']

        # States
        self.yaw_angle = config['yaw_angle_0']
        self.yaw_rate = config['yaw_rate_0']

        # Initialize
        self.positions, self.width, self.height, self.area, self.rows = self.calculate_hexagonal_positions()
        self.rotors = [Rotor(self.r, pos, self, self.q) for pos in self.positions]

    def update_all(self, delta_psi, delta_dpsi, selected_rotors=None, delta_pitch=None):
        self.yaw_angle += delta_psi
        self.yaw_rate += delta_dpsi
        
        if selected_rotors and delta_pitch:
            self.update_pitch_angles(selected_rotors, delta_pitch)

        self.update_all_rotors()

    def update_all_rotors(self):
        # Using list comprehension for efficiency
        [rotor.update_properties() for rotor in self.rotors]
    
    def update_pitch_angles(self, selected_rotors, delta_pitch):
        for index, dth in zip(selected_rotors, delta_pitch):
            self.rotors[index].pitch_angle += dth

    def calculate_hexagonal_positions(self):
        positions = []
        rows = []
        row_count = 0
        num_circles_placed = 0
        dx = 2 * self.r
        dy = self.r * np.sqrt(3)
        while num_circles_placed < self.n:
            row_length = self.x if row_count % 2 == 0 else self.x - 1
            if num_circles_placed + row_length > self.n:
                row_length = self.n - num_circles_placed
            offset = 0 if row_count % 2 == 0 else self.r
            row_positions = [(dx * i + offset + self.r - (self.x / 2 * dx), dy * row_count + self.r) for i in range(row_length)]
            positions.extend(row_positions)
            rows.append(row_positions)
            num_circles_placed += row_length
            row_count += 1
        width = max(pos[0] for pos in positions) + self.r
        height = max(pos[1] for pos in positions) + self.r
        area = width * height
        return positions, width, height, area, rows

    def calculate_ddpsi(self):
        total_moments = np.array([rotor.moments for rotor in self.rotors])
        total_thrust = sum(rotor.thrust for rotor in self.rotors)
        
        total_x_moment = np.sum(total_moments[:, 0])
        total_z_moment = np.sum(total_moments[:, 1])
        
        M_f = - np.sign(self.yaw_rate) * self.mu * (4.4 * total_x_moment * 1.356 + 
                self.W * 4.44822 * 2 * self.R * 0.3048 +
                2.2 * total_thrust * 4.44822 * 2 * self.R * 0.3048) / 2
        ddpsi = (total_z_moment + M_f ) / self.J
        return total_z_moment, M_f, ddpsi
        
    def select_half_array(self, half='left'):
        selected_indices = []
        for row in self.rows:
            row_indices = [self.positions.index(pos) for pos in row]
            if half == 'left':
                selected_indices.extend(row_indices[:len(row_indices) // 2])
            elif half == 'right':
                selected_indices.extend(row_indices[(len(row_indices) + 1) // 2:])
        return selected_indices
    
    def select_extremities(self, half='left'):
        left_extremities = []
        right_extremities = []
        for row in self.rows:
            if row:
                row_indices = [self.positions.index(pos) for pos in row]
                left_extremities.append(row_indices[0])
                right_extremities.append(row_indices[-1])
        
        if half == 'left':
            return left_extremities
        elif half == 'right':
            return right_extremities
    
    def plot_grid(self, selected_rotors=[]):
        """Plot circles with given positions and display thrust values."""
        fig, ax = plt.subplots()
        ax.set_xlim(-self.width, self.width)
        ax.set_ylim(0, self.height)
        ax.set_aspect('equal', 'box')
        
        for i, rotor in enumerate(self.rotors):
            x, y = rotor.position
            if i in selected_rotors:
                circle = plt.Circle((x, y), self.r / 1.05, facecolor='red', alpha=0.6)
            else:
                circle = plt.Circle((x, y), self.r / 1.05, facecolor='blue', alpha=0.6)
            ax.add_patch(circle)
            ax.text(x, y, f'T = {rotor.thrust:.2f}\n M = {rotor.moments[1]:.2f}', fontsize=6, ha='center')

        ax.set_title('Honeycomb')
        ax.set_xlabel("Width [m]")
        ax.set_ylabel("Height [m]")
        plt.show()