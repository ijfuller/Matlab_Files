clear all
clc

Receiver_Height=12
N_Tubes=10
D_tube=40/1000
Spacing_Factor=1;

One_Pass_Height=(N_Tubes*D_tube)+((N_Tubes-1)*(Spacing_Factor*D_tube))
N_Passes=floor(Receiver_Height/One_Pass_Height)

Q_Panel_Solar_Incident=10e6/4
eta=0.95;

Q_absorbed=Q_Panel_Solar_Incident*eta
Delta_T=432.5;
cp=1.522; %kJ/kg-K
rho=1813; %kg/m^3

m_dot=(Q_absorbed/1000)/(cp*Delta_T)
V_dot=m_dot/rho

m_dot_per_tube=m_dot/N_Tubes
V_dot_per_tube=V_dot/N_Tubes