function [h_bar] = S_and_K_Function(Avg_Temp, H, T_infinity)

g=9.81; %Gravity
L_c=H; %Characteristic length
T_avg=Avg_Temp;


alpha=4.889368666005050E-17*T_avg^4-1.794474675669010E-13*T_avg^3+2.98380685351056E-10*T_avg^2-1.00075949101992E-9*T_avg-3.44798307584351e-7;
nu=1.03450643178104E-17*T_avg^4-4.85019754418772E-14*T_avg^3+1.3580075963433E-10*T_avg^2+2.27985665430374E-8*T_avg-2.0313337298359E-6;
k=-1.24607229972985E-16*T_avg^4+5.01096786429384E-12*T_avg^3-2.940474355754410E-8*T_avg^2+9.05978900277077E-5*T_avg+9.82003734668099E-4;
beta=1/T_avg;

Ra=((g*beta)/(nu*alpha))*(T_avg-T_infinity)*L_c^3;
Gr_L=(g*beta*(T_avg-T_infinity)*L_c^3)/(nu^2);
Nu_L=0.088*Gr_L^(1/3)*(T_avg/T_infinity)^0.18;
h_bar=(Nu_L*k)/L_c;

if Gr_L<=10^5|| 10^12<=Gr_L
    %display('WARNING! Operation conditions out of correlation range!');
end
end