import power_curve_2 as P
import numpy as np
import xarray as x
import matplotlib.pyplot as plt
#print(P.P_array)

layoutreversed = np.array([5,6,5,6,5,6])
currentlayout = np.array([6,5,6,5,6,5])

height_array = np.array([45.8845, 82.05752, 118.2305, 154.4035, 190.5766, 226.7496])
# power_curve = P.P_array/33

power_curve = 0.5 * P.rho * P.CP * P.U_array ** 3 * np.pi * P.radius ** 2 
power_curve[power_curve > 30/33*10**6] = 30/33*10**6
power_curve[P.U_array < P.cut_in] = 0
power_curve[P.U_array > P.cut_off] = 0

# speeds_height_1 = np.zeros(8760)
# speeds_height_2 = np.zeros(8760)
# speeds_height_3 = np.zeros(8760)
# speeds_height_4 = np.zeros(8760)
# speeds_height_5 = np.zeros(8760)
# speeds_height_6 = np.zeros(8760)

# currentlayout[i]*

from scipy.interpolate import interp1d
    

power_function = interp1d(P.U_array, power_curve, kind='linear', fill_value="extrapolate") #for 0.91 MW generator
# Interpolate missing values
U_interp = np.linspace(min(P.U_array), max(P.U_array), num=50)  # Interpolated points
power_interp = power_function(U_interp)
# plt.scatter(P.U_array, power_curve, label='Original data with NaNs')
# plt.plot(U_interp, power_interp, label='Interpolated data', color='red')
# plt.legend()
# plt.xlabel('pitch')
# plt.ylabel('CT')
# plt.title('Interpolation of CT')
# plt.show()

# print(power_function(13.5))

# energy_1 = power_function(speeds_height_1)
# print(np.sum(energy_1))

import pickle

# Load the arrays from the file
with open('arrays.pkl', 'rb') as f:
    arrays = pickle.load(f)

# Access the arrays
array1, array2, array3, array4, array5, array6 = arrays

print("Arrays loaded successfully!")

# Calculate the average of each 2x3 matrix
speeds_height_1 = np.mean(array1, axis=(1, 2))
speeds_height_2 = np.mean(array2, axis=(1, 2))
speeds_height_3 = np.mean(array3, axis=(1, 2))
speeds_height_4 = np.mean(array4, axis=(1, 2))
speeds_height_5 = np.mean(array5, axis=(1, 2))
speeds_height_6 = np.mean(array6, axis=(1, 2))

t = np.arange(0, 8760, 1)

# print(speeds_height_1.shape)  # Should print (8760,)
# print(speeds_height_1)
# print(array1)
# plt.plot(t, speeds_height_4)
# plt.show()

power_curve_084 = 0.5 * P.rho * P.CP * P.U_array ** 3 * np.pi * P.radius ** 2 
power_curve_084[power_curve_084 > 0.84*10**6] = 0.84*10**6
power_curve_084[P.U_array < P.cut_in] = 0
power_curve_084[P.U_array > P.cut_off] = 0
power_function_2 = interp1d(P.U_array, power_curve_084, kind='linear', fill_value="extrapolate")  #0.84

#0.74, 0.82, 0,87, 0.92, 0.95, 0.98
power_curve_097 = 0.5 * P.rho * P.CP * P.U_array ** 3 * np.pi * P.radius ** 2
power_curve_097[power_curve_097 > 0.97*10**6] = 0.97*10**6
power_curve_097[P.U_array > P.cut_off] = 0
power_curve_097[P.U_array < P.cut_in] = 0
power_function_3 = interp1d(P.U_array, power_curve_097, kind='linear', fill_value="extrapolate")  #0.97

energy_1 = np.sum(power_function(speeds_height_1)*currentlayout[0]) + np.sum(power_function(speeds_height_2)*currentlayout[1])+ np.sum(power_function(speeds_height_3)*currentlayout[2])+ np.sum(power_function(speeds_height_4)*currentlayout[3])+ np.sum(power_function(speeds_height_5)*currentlayout[4])+ np.sum(power_function(speeds_height_6)*currentlayout[5]) # one size fits all config
installedpower1 = np.max(power_function(U_interp)*33)
energy_2 = np.sum(power_function_2(speeds_height_1)*currentlayout[0]) + np.sum(power_function_2(speeds_height_2)*currentlayout[1])+ np.sum(power_function_2(speeds_height_3)*currentlayout[2])+ np.sum(power_function_3(speeds_height_4)*currentlayout[3])+ np.sum(power_function_3(speeds_height_5)*currentlayout[4])+ np.sum(power_function_3(speeds_height_6)*currentlayout[5])  #3,3 config
installedpower2 = np.max(power_function_2(U_interp)*(6+5+6))+np.max(power_function_3(U_interp)*(5+6+5))



power_curve_076 = 0.5 * P.rho * P.CP * P.U_array ** 3 * np.pi * P.radius ** 2 
power_curve_076[power_curve_076 > 0.76*10**6] = 0.76*10**6
power_curve_076[P.U_array > P.cut_off] = 0
power_curve_076[P.U_array < P.cut_in] = 0
power_function_4 = interp1d(P.U_array, power_curve_076, kind='linear', fill_value="extrapolate")  #0.76

power_curve_086 = 0.5 * P.rho * P.CP * P.U_array ** 3 * np.pi * P.radius ** 2 
power_curve_086[power_curve_086 > 0.86*10**6] = 0.86*10**6
power_curve_086[P.U_array > P.cut_off] = 0
power_curve_086[P.U_array < P.cut_in] = 0
power_function_5 = interp1d(P.U_array, power_curve_086, kind='linear', fill_value="extrapolate")  #0.86
plt.plot(P.U_array, power_curve_084)
plt.show()
plt.plot(P.U_array, power_curve_097)
plt.show()
plt.plot(P.U_array, power_curve_076)
plt.show()
plt.plot(P.U_array, power_curve_086)
plt.show()
energy_3 = np.sum(power_function_4(speeds_height_1)*currentlayout[0]) + np.sum(power_function_5(speeds_height_2)*currentlayout[1])+ np.sum(power_function_5(speeds_height_3)*currentlayout[2])+ np.sum(power_function_3(speeds_height_4)*currentlayout[3])+ np.sum(power_function_3(speeds_height_5)*currentlayout[4])+ np.sum(power_function_3(speeds_height_6)*currentlayout[5])  #1,2,3 config
installedpower3 = np.max(power_function_4(U_interp)*6)+np.max(power_function_5(U_interp)*(5+6))+np.max(power_function_3(U_interp)*(5+6+5))


print(f'{energy_1=},{energy_2=},{energy_3=}') 

power_1 = energy_1/8760
power_2=energy_2/8760
power_3 = energy_3/8760
print(f'{power_1=},{power_2=},{power_3=}') 

plt.plot(P.U_array, power_curve)
plt.show()
print(f'{installedpower1=},{installedpower2=},{installedpower3=}') 

normalisedpower1 = power_1/installedpower1
normalisedpower2 = power_2/installedpower2
normalisedpower3 = power_3/installedpower3
print(f'{normalisedpower1=},{normalisedpower2=},{normalisedpower3=}') 

#data analysis now for frederico 
increase21 = normalisedpower2/normalisedpower1-1
increase31 = normalisedpower3/normalisedpower1-1
print(f'{increase21=}, {increase31=}')

averagewindspeed1 = np.max([arrays])

print(f'{averagewindspeed1=}')