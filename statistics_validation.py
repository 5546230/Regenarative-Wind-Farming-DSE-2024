import matplotlib.pyplot as plt
import numpy as np

def statistics():
    references = {}

    'IEA 22 MW'
    dict = {}
    dict['hub_height'] = 170
    dict['tower_mass'] = 1574
    dict['monopile_mass'] = 2097
    dict['max_t'] = 0.091
    dict['min_t'] = 0.038
    dict['max_D'] = 10
    dict['min_D'] = 6
    dict['h_above_mud'] = 34+170
    dict['m_top'] = (821.2+120+82*3)*1000
    dict['P_rated'] = 22e6
    dict['v_rated'] = 11
    dict['T_rated'] = 2793e3

    references['IEA22MW'] = dict


    '5 MW reference'
    dict = {}
    dict['hub_height'] = 90
    dict['tower_mass'] = 347460/1000
    dict['monopile_mass'] = 0
    dict['max_t'] = .027
    dict['min_t'] = .019
    dict['max_D'] = 6
    dict['min_D'] = 3.87
    dict['h_above_mud'] = 87.6
    dict['P_rated'] = 18.45e6
    dict['v_rated'] = 11.4
    dict['T_rated'] = 800e3

    dict['m_top'] = 110e3 + 240e3
    references['5MWref'] = dict

    'IEA 15 MW'
    dict = {}
    dict['hub_height'] = 150
    dict['tower_mass'] = 860
    dict['monopile_mass'] = 1318
    dict['max_t'] = .055
    dict['min_t'] =.02
    dict['max_D'] = 10
    dict['min_D'] = 6.5
    dict['h_above_mud'] = 150+30
    dict['P_rated'] = 31.77#15e6
    dict['v_rated'] = 10.59
    dict['m_top'] = 1017e3
    dict['T_rated'] = 2.8e6
    references['IEA15MW'] = dict

    'DTU 8 MW'
    dict = {}
    dict['hub_height'] = 110
    dict['tower_mass'] = 558e3
    dict['monopile_mass'] = 0
    dict['max_t'] = .036
    dict['min_t'] = .022
    dict['max_D'] = 7.7
    dict['min_D'] = 5
    dict['h_above_mud'] = 110
    dict['P_rated'] = 34.2875e6
    dict['v_rated'] = 12.5
    dict['m_top'] = 90e3 + 285e3+35e3*3
    dict['T_rated'] = 2743e3
    references['DTU8MW'] = dict

    #print(references.keys())
    fig, axs = plt.subplots(1, 2, figsize=(10, 5),  layout='constrained')

    plus1 = []
    minus1 = []
    plus2 = []
    minus2 = []
    for key in references.keys():

        turbine = references[key]
        axs[0].plot(turbine['h_above_mud']*np.ones(2), [turbine['max_D'], turbine['min_D']], label=str(key))
        axs[1].plot(turbine['h_above_mud'] * np.ones(2), [turbine['max_t'], turbine['min_t']], label=str(key))

        #single_tower = SingleTower(length=turbine['h_above_mud'], material=Steel(), M_truss=turbine['m_top'], sum_M_RNA=0, F_T=turbine['T_rated'], )
        #R_single, t_single, _ = single_tower.sizing_analysis(verbose=False)
        #D = 2 * R_single
        #t = t_single
        #axs[0].plot(turbine['h_above_mud'], D, marker='x')
        #axs[1].plot(turbine['h_above_mud'], t, marker='x')
        #plus1.append(abs(turbine['max_D'] - D)/D)
        #minus1.append(abs(turbine['min_D'] - D)/D)

        #plus2.append(abs(turbine['max_t'] - t) / t)
        #minus2.append(abs(turbine['min_t'] - t) / t)

    #print(np.average(plus1)*100, np.average(minus1)*100)
    #print(np.average(plus2) * 100, np.average(minus2) * 100)
    axs[0].set_xlabel('cylinder height')
    axs[1].set_xlabel('cylinder height')

    axs[0].set_ylabel('Diameter [m]')
    axs[1].set_ylabel('thickness [m]')
    axs[0].legend()
    axs[1].legend()
    #axs[0].set_title(f'uncertainty: +{np.average(plus1)*100:.2f}%, -{np.average(minus1)*100:.2f}%')
    #axs[1].set_title(f'uncertainty: +{np.average(plus2) * 100:.2f}%, -{np.average(minus2) * 100:.2f}%')
    plt.show()

    uncertainty_range_D = [np.average(plus1), np.average(minus1)]
    uncertainty_range_t = [np.average(plus2) , np.average(minus2)]
    return uncertainty_range_D, uncertainty_range_t
