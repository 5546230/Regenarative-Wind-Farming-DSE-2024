import numpy as np
import generatorselection as g

number_turbines = np.array([15,14,14,14,14,13,13,13,13,12,12,12,12,11,11,11,11,10,10,10,10,9,9,9,6])
assumed_regeneration = 0.912
assumed_regeneration_speed = assumed_regeneration**(1/3)

unit_cf = 0.59 
availability = 0.9538
# reg2 = 

# print(np.sum(number_turbines))
# energy = 1*assumed_regeneration**(number_turbines-1)
# print(energy)

total_energy = 0

for i in number_turbines:
    for ix in range(i):
        total_energy += assumed_regeneration**(ix+1-1)

print(total_energy)
energy_produced = total_energy*30
print(f'{energy_produced=}')
energy_percentage = energy_produced/30/288
print(f'{energy_percentage=}')
wind_farm_cf = energy_percentage*unit_cf*availability
print(f'{wind_farm_cf=}')
print(wind_farm_cf/unit_cf/availability)


speeds = [g.speeds_height_1, g.speeds_height_15, g.speeds_height_2, g.speeds_height_25, g.speeds_height_3, g.speeds_height_35, g.speeds_height_4, g.speeds_height_45, g.speeds_height_5, g.speeds_height_55, g.speeds_height_6]
data = np.array(speeds)
average_speeds = np.mean(data, axis=0)
print(np.shape(average_speeds))
print(g.power_function(average_speeds[0])*34)

# energy_final = np.sum(power_function(speeds_height_1)*currentlayout2[0]) + np.sum(power_function(speeds_height_15)*currentlayout2[1])+ np.sum(power_function(speeds_height_2)*currentlayout2[2])+ np.sum(power_function(speeds_height_25)*currentlayout2[3])+ np.sum(power_function(speeds_height_3)*currentlayout2[4])+ np.sum(power_function(speeds_height_35)*currentlayout2[5])+np.sum(power_function(speeds_height_4)*currentlayout2[6])+np.sum(power_function(speeds_height_45)*currentlayout2[7])+np.sum(power_function(speeds_height_5)*currentlayout2[8])+np.sum(power_function(speeds_height_55)*currentlayout2[9])+np.sum(power_function(speeds_height_6)*currentlayout2[10])


rot_per_row = np.array([4,2,4,2,4,2,4,2,4,2,4])
print(speeds[0].shape)

g.power_function()

regeneration_factors = [assumed_regeneration_speed**(i) for i in range(16)]


"""

for x in range(len(average_speeds)):
    for i in number_turbines:
        for ix in range(i):
            current_speed_1 = g.speeds_height_1[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_15 = g.speeds_height_15[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_2 = g.speeds_height_2[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_25 = g.speeds_height_25[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_3 = g.speeds_height_3[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_35 = g.speeds_height_35[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_4 = g.speeds_height_4[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_45 = g.speeds_height_45[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_5 = g.speeds_height_5[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_55 = g.speeds_height_55[x]*assumed_regeneration_speed**(ix+1-1)
            current_speed_6 = g.speeds_height_6[x]*assumed_regeneration_speed**(ix+1-1)
            total_energy += g.power_function(current_speed_1)*g.currentlayout2[0] + g.power_function(current_speed_15)*g.currentlayout2[1]+ g.power_function(current_speed_2)*g.currentlayout2[2]+ g.power_function(current_speed_25)*g.currentlayout2[3]+ g.power_function(current_speed_3)*g.currentlayout2[4]+ g.power_function(current_speed_35)*g.currentlayout2[5]+g.power_function(current_speed_4)*g.currentlayout2[6]+g.power_function(current_speed_45)*g.currentlayout2[7]+g.power_function(current_speed_5)*g.currentlayout2[8]+g.power_function(current_speed_55)*g.currentlayout2[9]+g.power_function(current_speed_6)*g.currentlayout2[10]


print('capacity factor is: ', total_energy/(30e06*288*8760)*availability)
print('AEP is:', total_energy/(30e06*288*8760)/unit_cf)
"""