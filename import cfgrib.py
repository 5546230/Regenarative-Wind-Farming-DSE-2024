import cfgrib
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def read_grib_file_cfgrib(filename):
    # Open the GRIB file as an xarray dataset
    ds = xr.open_dataset(filename, engine='cfgrib')
    
    # Print the dataset information
    print(ds)

    # Extracting u10 and v10 data
    u10 = ds['u10'].values
    v10 = ds['v10'].values
    
    # Calculate wind speed
    V = np.sqrt(u10**2 + v10**2)
    
    # Get time information
    time = ds['time'].values
    
    return V, u10, v10, time

# Load file
filename = 'C:/Users/Marcos/Desktop/TU Delft/BACHELOR - YEAR 3/DSE/WindDataDSEFinal.grib'
V, u10, v10, time = read_grib_file_cfgrib(filename)

# Define heights from 10m to 400m
heights = np.arange(10, 401, 1)  # heights in meters
z_r = 10  # reference height in meters
alpha = 0.06  # wind shear exponent

# Initialize a matrix to store V at different heights for all locations and times
matrix = np.zeros((V.shape[0], V.shape[1], V.shape[2], len(heights)))

# Adjust V for different heights
for t in range(V.shape[0]):
    for lat in range(V.shape[1]):
        for lon in range(V.shape[2]):
            for j, h in enumerate(heights):
                matrix[t, lat, lon, j] = V[t, lat, lon] * (h / z_r) ** alpha

# Calculate the average wind speed at each height over the entire dataset
V_height_avg = np.mean(matrix, axis=(0, 1, 2))

# Plotting the average wind speed at different heights
plt.figure(figsize=(12, 6))
plt.plot(V_height_avg, heights)
plt.title('Average Wind Speed at Different Heights')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Height (m)')
plt.grid(True)
plt.show()

# Adjust V for the specified height (150m)
V_150m = V * (150 / z_r) ** alpha

# Convert the wind speed data and time into a Pandas DataFrame
df = pd.DataFrame({
    'wind_speed': V_150m.flatten(),
    'time': pd.to_datetime(np.repeat(time, V_150m.shape[1] * V_150m.shape[2]))
})

# Resample to daily average
daily_avg = df.resample('D', on='time').mean()

# Resample to monthly average
monthly_avg = df.resample('M', on='time').mean()

# Calculate the mean wind speed
mean_speed = df['wind_speed'].mean()

# Plotting the daily, monthly averages and mean wind speed
plt.figure(figsize=(12, 6))
plt.plot(daily_avg.index, daily_avg['wind_speed'], label='Daily average')
plt.plot(monthly_avg.index, monthly_avg['wind_speed'], 'r-o', label='Monthly average')
plt.axhline(y=mean_speed, color='k', linestyle='--', label=f'Mean speed: {mean_speed:.2f} m/s')

# Formatting the x-axis to show month names
plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%b'))
plt.gca().xaxis.set_major_locator(plt.matplotlib.dates.MonthLocator())

plt.title('Wind speed measurements in 2023 at 150m')
plt.xlabel('Date')
plt.ylabel('Wind speed (m/s)')
plt.legend()
plt.grid(True)
plt.show()

# Creating bins for the wind rose
angle = (np.degrees(np.arctan2(v10, u10)) + 180) % 360
num_bins = 36  # 10-degree bins
bins = np.linspace(0, 360, num_bins + 1)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Counting the number of occurrences in each bin
hist, _ = np.histogram(angle.flatten(), bins)

# Plotting the wind rose
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))

# Convert angles from degrees to radians for the polar plot
theta = np.deg2rad(bin_centers)
width = np.deg2rad(bins[1] - bins[0])

bars = ax.bar(theta, hist, width=width, edgecolor='black')

# Add labels and title
ax.set_theta_zero_location('N')  # Set 0 degrees to the top (North)
ax.set_theta_direction(-1)  # Clockwise
ax.set_title('Wind Frequency Rose (Wind Coming From)')

plt.show()
