clear all
clc

g=9.81;
L_c=1;
T_infinity=293.15;
T_wall=800;
T_avg=(T_wall+T_infinity)/2


alpha=4.889368666005050E-17*T_avg^4-1.794474675669010E-13*T_avg^3+2.98380685351056E-10*T_avg^2-1.00075949101992E-9*T_avg-3.44798307584351e-7;
nu=1.03450643178104E-17*T_avg^4-4.85019754418772E-14*T_avg^3+1.3580075963433E-10*T_avg^2+2.27985665430374E-8*T_avg-2.0313337298359E-6;
k=-1.24607229972985E-16*T_avg^4+5.01096786429384E-12*T_avg^3-2.940474355754410E-8*T_avg^2+9.05978900277077E-5*T_avg+9.82003734668099E-4;
beta=1/T_avg;

Ra=((g*beta)/(nu*alpha))*(T_wall-T_infinity)*L_c^3
Gr_L=(g*beta*(T_wall-T_infinity)*L_c^3)/(nu^2)
Nu_L=0.088*Gr_L^(1/3)*(T_wall/T_infinity)^0.18
h_conv_0=(Nu_L*k)/L_c

display('-----------------------')
%%%mods with lip%%%
%%%cube model for simplicity L=W=H=2
H=2
W=2

Top_Lip=0
Bottom_Lip=0

if inclination_angle <=30an
    n=0.63;
else
    n=0.8;
end

Modifier=(A_1/A_2)*(A_3/A_1)^n
