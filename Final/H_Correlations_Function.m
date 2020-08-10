function [h_bar_conv] = H_Correlations_Function(Receiver_Height,T_h,A,T_infinity,h_Type)
    gravity=9.81;
    Lc=Receiver_Height;
    A_h=A;
    if length(A_h)==7
        A_h(7)=[];
    elseif length(A_h)==8
        A_h(8)=[];
    else
        A_h(9)=[];
    end
    A_Total=sum(A_h);

    T_avg=0;
    for i=1:length(A_h)
        T_avg=T_avg+T_h(i)*(A_h(i)/A_Total);
    end
    %%%Siebers and Kraabel"
    if h_Type==2
        T_film=(T_avg+T_infinity)/2;
        beta=1/T_infinity;
        nu=1.03450643178104E-17*T_infinity^4-4.85019754418772E-14*T_infinity^3+1.3580075963433E-10*T_infinity^2+2.27985665430374E-8*T_infinity-2.0313337298359E-6;
        k=-1.24607229972985E-16*T_infinity^4+5.01096786429384E-12*T_infinity^3-2.940474355754410E-8*T_infinity^2+9.05978900277077E-5*T_infinity+9.82003734668099E-4;
        Gr=(gravity*beta*(T_avg-T_infinity)*Lc^3)/nu^2;
        Nuss=0.088*Gr^(1/3)*(T_avg/T_infinity)^0.18;
        h_bar_conv=(Nuss*k)/Lc;
    end

end

