import power_curve_2 as P
import numpy as np
import xarray as x
import matplotlib.pyplot as plt
#print(P.P_array)

layoutreversed = np.array([5,6,5,6,5,6])
currentlayout = np.array([6,5,6,5,6,5])
currentlayout2 = np.array([4,2,4,2,4,2,4,2,4,2,4])

height_array = np.array([45.8845, 82.05752, 118.2305, 154.4035, 190.5766, 226.7496])
# power_curve = P.P_array/33

power_curve = 0.5 * P.rho * P.CP * P.U_array ** 3 * np.pi * P.radius ** 2 
power_curve[power_curve > 30/34*10**6] = 30/34*10**6
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
array1, array15,array2,array25, array3,array35, array4,array45, array5,array55, array6 = arrays

print("Arrays loaded successfully!")

# Calculate the average of each 2x3 matrix
speeds_height_1 = np.mean(array1, axis=(1, 2))
speeds_height_15 = np.mean(array15, axis=(1, 2))
speeds_height_2 = np.mean(array2, axis=(1, 2))
speeds_height_25 = np.mean(array25, axis=(1, 2))
speeds_height_3 = np.mean(array3, axis=(1, 2))
speeds_height_35 = np.mean(array35, axis=(1, 2))
speeds_height_4 = np.mean(array4, axis=(1, 2))
speeds_height_45 = np.mean(array45, axis=(1, 2))
speeds_height_5 = np.mean(array5, axis=(1, 2))
speeds_height_55 = np.mean(array55, axis=(1, 2))
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

energy_final = np.sum(power_function(speeds_height_1)*currentlayout2[0]) + np.sum(power_function(speeds_height_15)*currentlayout2[1])+ np.sum(power_function(speeds_height_2)*currentlayout2[2])+ np.sum(power_function(speeds_height_25)*currentlayout2[3])+ np.sum(power_function(speeds_height_3)*currentlayout2[4])+ np.sum(power_function(speeds_height_35)*currentlayout2[5])+np.sum(power_function(speeds_height_4)*currentlayout2[6])+np.sum(power_function(speeds_height_45)*currentlayout2[7])+np.sum(power_function(speeds_height_5)*currentlayout2[8])+np.sum(power_function(speeds_height_55)*currentlayout2[9])+np.sum(power_function(speeds_height_6)*currentlayout2[10]) # one size fits all config
print(f'{energy_final=}')
energy_installed = 30*10**6*8760
capacity_Factor = energy_final/energy_installed
print(f'{capacity_Factor=}')
print(np.max(power_function(U_interp)))

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
# plt.plot(P.U_array, power_curve_084)
# plt.show()
# plt.plot(P.U_array, power_curve_097)
# plt.show()
# plt.plot(P.U_array, power_curve_076)
# plt.show()
# plt.plot(P.U_array, power_curve_086)
# plt.show()
energy_3 = np.sum(power_function_4(speeds_height_1)*currentlayout[0]) + np.sum(power_function_5(speeds_height_2)*currentlayout[1])+ np.sum(power_function_5(speeds_height_3)*currentlayout[2])+ np.sum(power_function_3(speeds_height_4)*currentlayout[3])+ np.sum(power_function_3(speeds_height_5)*currentlayout[4])+ np.sum(power_function_3(speeds_height_6)*currentlayout[5])  #1,2,3 config
installedpower3 = np.max(power_function_4(U_interp)*6)+np.max(power_function_5(U_interp)*(5+6))+np.max(power_function_3(U_interp)*(5+6+5))


# print(f'{energy_1=},{energy_2=},{energy_3=}') 

power_1 = energy_1/8760
power_2=energy_2/8760
power_3 = energy_3/8760
# print(f'{power_1=},{power_2=},{power_3=}') 

# plt.plot(P.U_array, power_curve)
# plt.show()
# print(f'{installedpower1=},{installedpower2=},{installedpower3=}') 

normalisedpower1 = power_1/installedpower1
normalisedpower2 = power_2/installedpower2
normalisedpower3 = power_3/installedpower3
# print(f'{normalisedpower1=},{normalisedpower2=},{normalisedpower3=}') 

#data analysis now for frederico 
increase21 = normalisedpower2/normalisedpower1-1
increase31 = normalisedpower3/normalisedpower1-1
# print(f'{increase21=}, {increase31=}')

averagewindspeed1 = np.average([speeds_height_1])
averagewindspeed15 = np.average([speeds_height_15])
averagewindspeed2 = np.average([speeds_height_2])
averagewindspeed25 = np.average([speeds_height_25])
averagewindspeed3 = np.average([speeds_height_3])
averagewindspeed35 = np.average([speeds_height_35])
averagewindspeed4 = np.average([speeds_height_4])
averagewindspeed45 = np.average([speeds_height_45])
averagewindspeed5 = np.average([speeds_height_5])
averagewindspeed55 = np.average([speeds_height_55])
averagewindspeed6 = np.average([speeds_height_6])
AverageWindSpeed = np.average(arrays)
print(f'{AverageWindSpeed=}')
# print(f'{averagewindspeed6=}')
# arr = np.array([1, 3, 5, 7, 9])
count = np.sum(speeds_height_6 > 12)
count_percentage = count/8760
# print(f'{count_percentage=}')
# print(f'{averagewindspeed1=},{averagewindspeed15=},{averagewindspeed2=},{averagewindspeed25=},{averagewindspeed3=},{averagewindspeed35=},{averagewindspeed4=},{averagewindspeed45=},{averagewindspeed5=},{averagewindspeed55=},{averagewindspeed6=},')
# print(P.CT)
AREA = P.AREA/P.n_rotors
thrust1 = P.CT*0.5*P.rho*averagewindspeed1**2*AREA
thrust15 = P.CT*0.5*P.rho*averagewindspeed15**2*AREA
thrust2 = P.CT*0.5*P.rho*averagewindspeed2**2*AREA
thrust25 = P.CT*0.5*P.rho*averagewindspeed25**2*AREA
thrust3 = P.CT*0.5*P.rho*averagewindspeed3**2*AREA
thrust35 = P.CT*0.5*P.rho*averagewindspeed35**2*AREA
thrust4 = P.CT*0.5*P.rho*averagewindspeed4**2*AREA
thrust45 = P.CT*0.5*P.rho*averagewindspeed45**2*AREA
thrust5 = P.CT*0.5*P.rho*averagewindspeed5**2*AREA
thrust55 = P.CT*0.5*P.rho*averagewindspeed55**2*AREA
thrust6 = P.CT*0.5*P.rho*averagewindspeed6**2*AREA
# print(f'{thrust1=},{thrust15=},{thrust2=},{thrust25=},{thrust3=},{thrust35=},{thrust4=},{thrust45=},{thrust5=},{thrust55=},{thrust6=},')
# CT_1 = P.BEM.BEMspeedpotter(speeds_height_1)
# CT_15 = P.BEM.BEMspeedplotter(speeds_height_15)
# CT_2 = P.BEM.BEMspeedplotter(speeds_height_2)
# CT_25 = P.BEM.BEMspeedplotter(speeds_height_25)
# CT_3 = P.BEM.BEMspeedplotter(speeds_height_3)
# CT_35 = P.BEM.BEMspeedplotter(speeds_height_35)
# CT_4 = P.BEM.BEMspeedplotter(speeds_height_4)
# CT_45 = P.BEM.BEMspeedplotter(speeds_height_45)
# CT_5 = P.BEM.BEMspeedplotter(speeds_height_5)
# CT_55 = P.BEM.BEMspeedplotter(speeds_height_55)
# CT_6 = P.BEM.BEMspeedplotter(speeds_height_6)

# print(f'{speeds_height_1=},{speeds_height_15=},{speeds_height_2=},{speeds_height_25=},{speeds_height_3=},{speeds_height_35=},{speeds_height_4=},{speeds_height_45=},{speeds_height_5=},{speeds_height_55=},{speeds_height_6=},')