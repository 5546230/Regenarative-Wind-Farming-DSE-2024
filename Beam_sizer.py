#Beam sizer


W_truss = 4355e3*9.80665
W_RNA = 2602578.544*9.80665
L_AFC = 6.5e6
delta_max = 0.1
R_tower = 14.18/2
d = 27.5
P = (W_truss+W_RNA+L_AFC)
n_contact = 12
I_req = P/n_contact*((d-R_tower)**3)/(3*190e9*delta_max)
print(I_req)
rho = 7850


#=========== FIND I OF CHINESE BB7 BEAM ===============

A = 2700/rho
t_flange = (A-0.182)/(1.02-0.28)
print(t_flange)


w_flange = 0.51
h = 1.300
t_web = 0.14

I_cn = 2*((w_flange*(t_flange**3)/12)+(w_flange*t_flange*((h/2-t_flange/2)**2)))+(t_web*((h-2*t_flange)**3)/12)
print(I_cn)

delta_actual  = P/n_contact*((d)**3)/(3*190e9*I_cn)
sigma_actual = P*d/n_contact*0.65/I_cn
print(delta_actual, sigma_actual)



