import numpy as np
import generatorselection as g
import matplotlib.pyplot as plt

number_turbines = np.array([15,14,14,14,14,13,13,13,13,12,12,12,12,11,11,11,11,10,10,10,10,9,9,9,6])
regenerationsss = np.arange(0.2, 0.96, 0.05)
capacityfactorsss = []


for i in regenerationsss:
    assumed_regeneration = i #0.95 would be ideal but someone chose the wrong number of wings. 
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

    # g.power_function()

    # regeneration_factors = [assumed_regeneration_speed**(i) for i in range(16)]
    # regeneration_factorss = []
    # for i in (number_turbines):
    #     regeneration_factorss.append(regeneration_factors[:i])
    # # print()

    # print(regeneration_factors)
    speeds_heights = [
        g.speeds_height_1, g.speeds_height_15, g.speeds_height_2, g.speeds_height_25, 
        g.speeds_height_3, g.speeds_height_35, g.speeds_height_4, g.speeds_height_45, 
        g.speeds_height_5, g.speeds_height_55, g.speeds_height_6
    ]
    speeds_heights = np.array(speeds_heights)  # Shape: (11, 8760)

    # Initialize total_energy
    total_energy = 0

    # Create an array of regeneration factors
    regeneration_factors = np.power(assumed_regeneration_speed, np.arange(max(number_turbines)))
    # regeneration_factorss=[]
    # rege
    # for i in (number_turbines):
    #     #  print(regeneration_factors[:i])
    #     #  np.hstack((regeneration_factorss, regeneration_factors[:i]))
    #     regeneration_factorss.append(regeneration_factors[:i])
    #     for ix in range(15-i):
    #         regeneration_factorss.append(0)
    regeneration_factorss = np.zeros((15,len(number_turbines)))
    for i in range(len(number_turbines)):
        for ix in range(number_turbines[i]):
            regeneration_factorss[ix,i] = regeneration_factors[ix]




    # print(regeneration_factorss)
    # regeneration_factorss = np.array(regeneration_factorss)
    print(regeneration_factors)
    print(regeneration_factorss)

    # for x in range(len(average_speeds)):
    #     # Calculate the current speeds for all turbines
    #     current_speeds = speeds_heights[:, x].reshape(-1, 1) * regeneration_factors[:max(number_turbines)].reshape(1, -1)
    #     # print(regeneration_factors[:max(number_turbines)])

    #     # Calculate the power for each turbine speed
    #     power_matrix = g.power_function(current_speeds)

    #     # Multiply power by the corresponding layout factor for each speed height
    #     energy_contributions = power_matrix * g.currentlayout2[:len(speeds_heights)].reshape(-1, 1)

    #     # Sum the contributions to the total energy
    #     total_energy += np.sum(energy_contributions)

    # print(total_energy)
    # print('capacity factor is: ', total_energy/(30e06*288*8760)*availability)
    # print('AEP is:', total_energy/(30e06*288*8760)/unit_cf)

    # for x in range(len(average_speeds)):
    #     for i in number_turbines:
    #         for ix in range(i):
    #             current_speed = average_speeds[x]*assumed_regeneration_speed**(ix+1-1)
    #             total_energy += g.power_function(current_speed)*34 

    # print('capacity factor is: ', total_energy/(30e06*288*8760)*availability)
    # print('AEP is:', total_energy/(30e06*288*8760)/unit_cf)


    # for x in range(len(average_speeds)):
    #     for i in number_turbines:
    #         for ix in range(i):
    #             current_speed_1 = g.speeds_height_1[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_15 = g.speeds_height_15[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_2 = g.speeds_height_2[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_25 = g.speeds_height_25[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_3 = g.speeds_height_3[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_35 = g.speeds_height_35[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_4 = g.speeds_height_4[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_45 = g.speeds_height_45[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_5 = g.speeds_height_5[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_55 = g.speeds_height_55[x]*assumed_regeneration_speed**(ix+1-1)
    #             current_speed_6 = g.speeds_height_6[x]*assumed_regeneration_speed**(ix+1-1)
    #             total_energy += g.power_function(current_speed_1)*g.currentlayout2[0] + g.power_function(current_speed_15)*g.currentlayout2[1]+ g.power_function(current_speed_2)*g.currentlayout2[2]+ g.power_function(current_speed_25)*g.currentlayout2[3]+ g.power_function(current_speed_3)*g.currentlayout2[4]+ g.power_function(current_speed_35)*g.currentlayout2[5]+g.power_function(current_speed_4)*g.currentlayout2[6]+g.power_function(current_speed_45)*g.currentlayout2[7]+g.power_function(current_speed_5)*g.currentlayout2[8]+g.power_function(current_speed_55)*g.currentlayout2[9]+g.power_function(current_speed_6)*g.currentlayout2[10]

    # print('capacity factor is: ', total_energy/(30e06*288*8760)*availability)
    # print('AEP is:', total_energy/(30e06*288*8760)/unit_cf)
    energy = 0
    for y in range(len(number_turbines)):
        for x in range(15):
            energy+= np.sum(g.power_function(g.speeds_height_1*regeneration_factorss[x,y])*g.currentlayout2[0]) + np.sum(g.power_function(g.speeds_height_15*regeneration_factorss[x,y])*g.currentlayout2[1])+ np.sum(g.power_function(g.speeds_height_2*regeneration_factorss[x,y])*g.currentlayout2[2])+ np.sum(g.power_function(g.speeds_height_25*regeneration_factorss[x,y])*g.currentlayout2[3])+ np.sum(g.power_function(g.speeds_height_3*regeneration_factorss[x,y])*g.currentlayout2[4])+ np.sum(g.power_function(g.speeds_height_35*regeneration_factorss[x,y])*g.currentlayout2[5])+np.sum(g.power_function(g.speeds_height_4*regeneration_factorss[x,y])*g.currentlayout2[6])+np.sum(g.power_function(g.speeds_height_45*regeneration_factorss[x,y])*g.currentlayout2[7])+np.sum(g.power_function(g.speeds_height_5*regeneration_factorss[x,y])*g.currentlayout2[8])+np.sum(g.power_function(g.speeds_height_55*regeneration_factorss[x,y])*g.currentlayout2[9])+np.sum(g.power_function(g.speeds_height_6*regeneration_factorss[x,y])*g.currentlayout2[10]) # one size fits all config
    print('capacity factor is: ', energy/(30e06*288*8760)*availability)
    print('AEP efficiency is:', energy/(30e06*288*8760)/unit_cf)
    capacityfactorsss.append(energy/(30e06*288*8760)*availability)


plt.plot(regenerationsss,capacityfactorsss)
plt.grid()
plt.xlabel('Regeneration Factors')
plt.ylabel('Capacity Factors')
plt.show()
print(regenerationsss)
print(np.array(capacityfactorsss))
