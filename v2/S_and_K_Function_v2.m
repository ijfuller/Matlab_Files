function [h_bar] = S_and_K_Function_v2(Avg_Temp, H, T_infinity)

    g=9.81; %Gravity
    L_c=H; %Characteristic length
    T_avg=Avg_Temp;

    h_bar=0.81*(T_avg-T_infinity)^0.426;
end