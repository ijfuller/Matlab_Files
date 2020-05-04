clear all
clc

%Bring in Areas, Temps, Receiver Height, h_Type 

h_Type=1;

gravity=9.81;
Lc=12;
T=[700 700 700 700 700 700];
T_infinity=293.15;
A=[64.2908 64.2908 64.2908 64.2908 69.2965 69.2695];
A_Total=sum(A);

T_avg=0;
for i=1:length(A-1)
    T_avg=T_avg+T(i)*(A(i)/A_Total);
end

%%%Siebers and Kraabel"
if h_Type==1
    T_film=(T_avg+T_infinity)/2;
    beta=1/T_film;
    nu=1.03450643178104E-17*T_infinity^4-4.85019754418772E-14*T_infinity^3+1.3580075963433E-10*T_infinity^2+2.27985665430374E-8*T_infinity-2.0313337298359E-6;
    k=-1.24607229972985E-16*T_infinity^4+5.01096786429384E-12*T_infinity^3-2.940474355754410E-8*T_infinity^2+9.05978900277077E-5*T_infinity+9.82003734668099E-4;
    Gr=(gravity*beta*(T_avg-T_infinity)*Lc^3)/nu^2;
    Nuss=0.088*Gr^(1/3)*(T_avg/T_infinity)^0.18;
    h_bar_conv=(Nuss*k)/Lc
end


%%%Clausing 1987%%%
if h_Type==2
    T_film=(T_avg+T_infinity)/2
    beta=1/T_film;
    nu_f=1.03450643178104E-17*T_film^4-4.85019754418772E-14*T_film^3+1.3580075963433E-10*T_film^2+2.27985665430374E-8*T_film-2.0313337298359E-6;
    k_f=-1.24607229972985E-16*T_film^4+5.01096786429384E-12*T_film^3-2.940474355754410E-8*T_film^2+9.05978900277077E-5*T_film+9.82003734668099E-4;
    Pr_f=-3.58027612400237E-17*T_film^5+3.12318964774912E-13*T_film^4-9.62623685786222E-10*T_film^3+1.36413319370172E-6*T_film^2-8.4947006070093E-4*T_film+8.82699497542028E-1;
    Pr_infinity=-3.58027612400237E-17*T_infinity^5+3.12318964774912E-13*T_infinity^4-9.62623685786222E-10*T_infinity^3+1.36413319370172E-6*T_infinity^2-8.4947006070093E-4*T_infinity+8.82699497542028E-1;
    Gr=(gravity*beta*(T_avg-T_infinity)*Lc^3)/nu_f^2;
    Ra=Gr*Pr_f;
    Ra_infinity=Gr*Pr_infinity;
    
    
    
    g=0.108*Ra^(1/3)
    f=0.2524+0.9163*(T_avg/T_infinity)-0.1663*(T_avg/T_infinity)^2
    
end
