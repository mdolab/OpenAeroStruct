span = 58.7630524 # [m] baseline CRM
W0 = 0.5 * 2.5e6 # [N] (MTOW of B777 is 3e5 kg with fuel)
CT = 9.81 * 17.e-6 # [1/s] (9.81 N/kg * 17e-6 kg/N/s)
R = 14.3e6 # [m] maximum range
M = 0.84 # at cruise

alpha = 5. # [deg.]
rho = 0.38 # [kg/m^3] at 35,000 ft
a = 295.4 # [m/s] at 35,000 ft
mu = 0.0000144446 # [N*s/m^2] at 35,000 ft
nu = mu / rho
Re = 5.e6 # 5 million, used for wind-tunnel and CFD comparisons
# Re = 40.e6 # 40 million, used for actual flow conditions over CRM wing
v = a * M

CL0 = 0.2
CD0 = 0.015
