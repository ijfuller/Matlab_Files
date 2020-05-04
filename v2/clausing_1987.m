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
b_fac=(T_wall-T_avg)/(T_wall-T_infinity)
Ra_lam=3.8*10^8;
Ra_turb=1.6*10^9;

if Ra <= 3.8*10^8;
    Flow_Type=1;
    Flow='Laminar'
    g_fac=0.63*Ra^0.25;
    f_fac=1;
elseif 3.8*10^8 < Ra && Ra < 1.6*10^9;
    Flow_Type=2;
    Flow='Transitional'
    g_fac=0.63*Ra^0.25
    f_fac_turb=0.2524+0.9163*(T_wall/T_infinity)-0.1663*(T_wall/T_infinity)^2;
    f_fac=(f_fac_turb-1)*((Ra^(1/3)-Ra_lam^(1/3))/(Ra_turb^(1/3)-Ra_lam^(1/3)))+1;   
else 1.6*10^9 <= Ra;
    Flow_Type=3;
    Flow='Turbulent'
    g_fac=0.108*Ra^0.25
    f_fac=0.2524+0.9163*(T_wall/T_infinity)-0.1663*(T_wall/T_infinity)^2
end

Nu=b_fac*g_fac*f_fac
h=(Nu*k)/L_c
