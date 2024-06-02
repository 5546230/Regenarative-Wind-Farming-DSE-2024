import numpy as np
import matplotlib.pyplot as plt
from HoneyComb import HoneycombGrid

config = {
    "n": 33,
    "r": 21,
    "x": 6,
    "q": 1.225 * 0.5 *  13.5,  # Example dynamic pressure
    "W": 6600*1000*9.81,    # Example weight
    "J": 1e9,     # Example moment of inertia
    "R": 5,      # Example radius for moment calculation
    "mu": 0.006,  # Example friction coefficient
    "dt" : 0.01,
    "tf" : 100,
    "yaw_angle_0" : (-10 * np.pi / 180),
    "yaw_rate_0" : 0,
}

class Simulation:
    def __init__(self, config, Kp, Ki, Kd):
        self.grid = HoneycombGrid(config)
        self.dt = config["dt"]
        self.times = np.arange(0, config["tf"], self.dt)
        self.yaw_angles = []
        self.yaw_rates = []
        self.pitch_angles = []
        self.pid_controller = PIDController(Kp, Ki, Kd, self.dt)
        self.delta_pitch_update_interval = 2  # Update interval in seconds

    def run(self):
        dpsi = self.grid.yaw_rate
        left_half_rotors = self.grid.select_half_array('left')
        right_half_rotors = self.grid.select_half_array('right')
        """"
        dtheta = -10  # Example pitch angle change
        delta_pitch = [dtheta] * len(left_half_rotors)
        self.grid.update_pitch_angles(left_half_rotors, delta_pitch)
        """
        for t in self.times:
            M_z, M_f, ddpsi = self.grid.calculate_ddpsi()
            delta_dpsi = ddpsi * self.dt
            delta_psi = dpsi * self.dt
            dpsi += delta_dpsi  # Update dpsi with the new delta_dpsi
            self.grid.update_all(delta_psi, delta_dpsi)

            self.yaw_angles.append(self.grid.yaw_angle)
            self.yaw_rates.append(dpsi)
            self.pitch_angles.append(self.grid.rotors[left_half_rotors[0]].pitch_angle)

            # Every 2 seconds, update the delta pitch using PID controller
            if t % self.delta_pitch_update_interval == 0:
                pitch_adjustment = self.pid_controller.control(self.grid.yaw_angle)
                delta_pitch = [pitch_adjustment] * len(left_half_rotors)                
                self.grid.update_pitch_angles(left_half_rotors, delta_pitch)
                print(f"Time: {t}, Pitching Moment: {M_z}, Friction Moment: {M_f}")
                print(f"Time: {t}, Pitch Adjustment: {pitch_adjustment}, Yaw Angle: {self.grid.yaw_angle}")

        self.plot_results(self.times, self.yaw_angles, self.yaw_rates, self.pitch_angles)

    def plot_results(self, times, yaw_angles, yaw_rates, pitch_angles):
        # Convert yaw angles, rates, and pitch angles to degrees
        yaw_angles_deg = np.rad2deg(yaw_angles)
        yaw_rates_deg_per_s = np.rad2deg(yaw_rates)
        pitch_angles_deg = pitch_angles

        # Filter the data to include only the times when yaw_angles are between -90 and 90 degrees
        mask = (yaw_angles_deg >= -90) & (yaw_angles_deg <= 90)
        filtered_times = times[mask]
        filtered_yaw_rates = yaw_rates_deg_per_s[mask]

        plt.figure(figsize=(15, 5))

        # Plot yaw angles over time
        plt.subplot(1, 3, 1)
        plt.plot(times, yaw_angles_deg)
        plt.title('Yaw Angle Over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('$\psi$ [deg]')
        plt.ylim(bottom=-90, top=90)
        plt.grid(True)

        # Plot yaw rates over time, filtered by yaw angles within the range [-90, 90] degrees
        plt.subplot(1, 3, 2)
        plt.plot(filtered_times, filtered_yaw_rates)
        plt.title('Yaw Angular Velocity Over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('$\dot{\psi}$ [deg/s]')
        plt.ylim(bottom=-10, top=10)  # Ensure the same y-range as the yaw angle plot
        plt.grid(True)

        # Plot pitch angles over time
        plt.subplot(1, 3, 3)
        plt.plot(times, pitch_angles_deg)
        plt.title('Pitch Angle Over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('$\\theta$ [deg]')
        plt.grid(True)

        # Show the plots
        plt.tight_layout()
        plt.show()

    
class PIDController:
    def __init__(self, Kp, Ki, Kd, dt):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.dt = dt
        self.y_desired = 0
        self.integral_error = 0
        self.previous_error = 0

    def control(self, y_current):
        error = self.y_desired - y_current
        self.integral_error += error * self.dt
        derivative_error = (error - self.previous_error) / self.dt
        output = self.Kp * error + self.Ki * self.integral_error + self.Kd * derivative_error
        self.previous_error = error
        return output

    def set_desired_y(self, y_desired):
        self.y_desired = y_desired

Kp = -1
Ki = 1
Kd = -0.8
    
# Initialize and run the simulation
sim = Simulation(config, Kp, Ki, Kd)
sim.run()