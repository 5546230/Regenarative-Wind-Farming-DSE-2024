import cfgrib
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_grib_file_cfgrib(filename):
    # Open the GRIB file as an xarray dataset
    ds = xr.open_dataset(filename, engine='cfgrib')

    # Extracting u10 and v10 data
    u10 = ds['u10'].values
    v10 = ds['v10'].values

    # Calculate wind speed
    V = np.sqrt(u10**2 + v10**2)

    # Get time information
    time = ds['time'].values

    return V, u10, v10, time


# Load file
filename = 'WindDataDSEFinal.grib'
V, u10, v10, time = read_grib_file_cfgrib(filename)



# Creating bins for the wind rose
angle = (np.degrees(np.arctan2(v10, u10)) + 180) % 360
windspeed = np.array(V)
bins = np.linspace(0, 360, 37)
hist, _ = np.histogram(angle, bins=bins, density=True)
centers = (bins[:-1] + bins[1:]) / 2

# create a DataFrame
df = pd.DataFrame({'angle': centers.flatten(), 'frequency': hist.flatten()})
df.to_csv('angle.csv', index=False)
