#Script to find the rated wind speed
'''
#APPROACH 1 - Using graph (k and V_h_mean can be calculated from APPROACH 2 below)
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def weibull_pdf(v, k, c):
    """ Calculate Weibull probability density function """
    return (k / c) * (v / c)**(k - 1) * np.exp(-(v / c)**k)

def scale_parameter(v_m, k):
    """ Calculate the scale parameter 'c' for the Weibull distribution """
    return v_m * k**2.6674 / (0.184 + 0.816 * k**2.73855)

def p_annual(v_cut_in, v_rated, v_cut_out, k, v_m):
    """ Calculate the annual power produced """
    c = scale_parameter(v_m, k)
    
    # Integral from cut-in speed to rated speed
    first_integral = quad(lambda v: v**3 * weibull_pdf(v, k, c), v_cut_in, v_rated)[0]
    
    # Integral from rated speed to cut-out speed
    second_integral = v_rated**3 * quad(lambda v: weibull_pdf(v, k, c), v_rated, v_cut_out)[0]
    
    return first_integral + second_integral

# Parameters
V_cut_in = 3  # cut-in speed in m/s
V_cut_out = 30  # cut-out speed in m/s
k = 2.3805 # shape parameter of Weibull distribution
V_h_mean = 9.249 # mean wind speed in m/s

# Generate range of rated speeds
V_rated_range = np.linspace(V_cut_in, V_cut_out, 100)
annual_powers = [p_annual(V_cut_in, V_rated, V_cut_out, k, V_h_mean) for v_rated in V_rated_range]

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(V_rated_range, annual_powers, label='Annual Power Output', color='blue')
plt.axvline(x=V_cut_in, linestyle='--', color='green', label='Cut-in Speed')
plt.axvline(x=V_cut_out, linestyle='--', color='purple', label='Cut-out Speed')
plt.title('Annual Power Output vs Rated Speed')
plt.xlabel('Rated Speed (m/s)')
plt.ylabel('Annual Power Output (MWh)')
plt.legend()
plt.grid(True)
plt.show()
'''


#APPROACH 2 - Assuming that the rated wind speed = wind speed carrying maximum energy
import xarray as xr
import numpy as np

def read_grib_file_cfgrib(filename):
    # Open the GRIB file as an xarray dataset
    ds = xr.open_dataset(filename, engine='cfgrib')
    
    # Print the dataset information
    print(ds)

    # Extracting u10 and v10 data
    u10 = ds['u10'].values
    v10 = ds['v10'].values
    V = np.sqrt(u10**2 + v10**2)
    
    return u10, v10, V

# Example usage
filename = 'C:/Users/Marcos/Desktop/TU Delft/BACHELOR - YEAR 3/DSE/WindDataDSEFinal.grib'
u10, v10, V = read_grib_file_cfgrib(filename)

# Define the height at which we want to adjust the wind speed
height = 150  # height in meters
z_r = 10  # reference height in meters
alpha = 0.06  # wind shear exponent (unstable = 0.06, neutral = 0.10, stable = 0.27)

# Adjust V for the specified height for all locations and times
V_h = V * (height / z_r) ** alpha

# Calculate the mean and standard deviation of the wind speed at 200m over all locations and times
V_h_mean = np.mean(V_h)
V_h_std = np.std(V_h)

# Estimate the shape factor k
k = (V_h_std / V_h_mean) ** -1.086
print(f"Estimated shape factor (k): {k:.4f}")

#Estimate scale parameter c
c = V_h_mean * (0.568 + 0.433 / k)**(-1 / k)
print(f"Estimated scale parameter (c): {c:.4f}")

#c2 = V_200m_mean * k**2.6674 / (0.184 + 0.816*k**2.73855)
#print(f"Estimated scale parameter (c2): {c2:.4f}")

#Calculate the wind speed carrying maximum energy (V_maxE)
V_maxE = c * (1 + 2/k)**(1/k)

#Assumption: wind speed carrying maximum energy = rated wind speed (V_maxE = V_R)
V_R = V_maxE
print(f'V_R = {V_R}')

#Calculate most probable speed V_mp
#V_mp = c * (1 - 1/k)**(1/k)
#print(f'V_mp= {V_mp}')

#print(f'V_h_mean = {V_h_mean}')
#print(f'k = {k}')
