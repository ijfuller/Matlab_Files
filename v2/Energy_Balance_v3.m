clear all
clc

%%%Solar emissivity
e_act_sol=0.95; %Emissivity in short wave range for active surfaces
e_pass_sol=0.95;%Emissivity in short wave range for passive surfaces
epsilon_sol=[e_act_sol;e_act_sol;e_act_sol;e_act_sol;e_pass_sol;e_pass_sol;1]; %Vector of solar emissivities
rho_sol=1-epsilon_sol; %Solar reflectivity vector
%%%Thermal emissivity
e_act_therm=0.4; %Emissivity in long wave range for active surfaces
e_pass_therm=0.4; %Emissivity in long wave range for passive surfaces
epsilon_therm=[e_act_therm;e_act_therm;e_act_therm;e_act_therm;e_pass_therm;e_pass_therm;1];%Vector of thermal emissivities
rho_therm=1-epsilon_therm;%Vector of thermal reflectivity

UA_HTF=1000; %Overall heat transfer coeffecient of HTF fluid
T_infinity=273.15+20;
T_avg_HTF=705.65;

N_el_domain=7;
VF_Matrix=zeros(N_el_domain,N_el_domain); %Initialize view factor full matrix
sigma=5.67*10^-8;

A_1=5.36*12;%Panel_Area
A_2=69.36; %Floor Area
A_3=12*14; %Aperature Area

A(1)=A_1;
A(2)=A_1;
A(3)=A_1;
A(4)=A_1;
A(5)=A_2;
A(6)=A_2;
A(7)=A_3;

%Panel 1
VF_Matrix(1,1)=0; %Panel 1 to panel 1
VF_Matrix(1,2)=0.06138; %Panel 1 to panel 2
VF_Matrix(1,3)=0.08354; %Panel 1 to panel 3
VF_Matrix(1,4)=0.0918; %Panel 1 to panel 4
VF_Matrix(1,7)=0.4822; %Panel 1 to opening
VF_Matrix(1,5)=(1-sum(VF_Matrix(1,:)))/2; %Panel 1 to floor
VF_Matrix(1,6)=VF_Matrix(1,5); %Panel 1 to roof
%Panel 2
VF_Matrix(2,1)=0.06138; %Panel 2 to panel 1
VF_Matrix(2,2)=0; %Panel 2 to panel 2
VF_Matrix(2,3)=0.06138; %Panel 2 to panel 3
VF_Matrix(2,4)=0.08354; %Panel 2 to panel 4
VF_Matrix(2,7)=0.4717; %Panel 2 to opening
VF_Matrix(2,5)=(1-sum(VF_Matrix(2,:)))/2; %Panel 2 to floor
VF_Matrix(2,6)=VF_Matrix(2,5); %Panel 2 to roof
%Panel 3
VF_Matrix(3,1)=0.08354; %Panel 3 to panel 1
VF_Matrix(3,2)=0.06138; %Panel 3 to panel 2
VF_Matrix(3,3)=0; %Panel 3 to panel 3
VF_Matrix(3,4)=0.06138; %Panel 3 to panel 4
VF_Matrix(3,7)=0.4717; %Panel 3 to opening
VF_Matrix(3,5)=(1-sum(VF_Matrix(3,:)))/2; %Panel 3 to floor
VF_Matrix(3,6)=VF_Matrix(3,5); %Panel 3 to roof
%Panel 4
VF_Matrix(4,1)=0.0918; %Panel 4 to panel 1
VF_Matrix(4,2)=0.08354; %Panel 4 to panel 2
VF_Matrix(4,3)=0.06138; %Panel 4 to panel 3
VF_Matrix(4,4)=0; %Panel 4 to panel 4
VF_Matrix(4,7)=0.4822; %Panel 4 to opening
VF_Matrix(4,5)=(1-sum(VF_Matrix(4,:)))/2; %Panel 4 to floor
VF_Matrix(4,6)=VF_Matrix(4,5); %Panel 4 to roof
%Floor
VF_Matrix(5,1)=(VF_Matrix(1,5)*A_1)/A_2; %Floor to panel 1
VF_Matrix(5,2)=(VF_Matrix(2,5)*A_1)/A_2; %Floor to panel 2
VF_Matrix(5,3)=(VF_Matrix(3,5)*A_1)/A_2; %Floor to panel 3
VF_Matrix(5,4)=(VF_Matrix(4,5)*A_1)/A_2; %Floor to panel 4
VF_Matrix(5,7)=(0.1350*A_3)/A_2; %Floor to opening
VF_Matrix(5,5)=0; %Floor to floor
VF_Matrix(5,6)=1-sum(VF_Matrix(5,:)); %floor to roof
%Roof
VF_Matrix(6,:)=VF_Matrix(5,:);
VF_Matrix(6,5)=VF_Matrix(6,6);
VF_Matrix(6,6)=0;
%Aperature
VF_Matrix(7,1)=(VF_Matrix(1,7)*A_1)/A_3; %%Aperature to panel 1
VF_Matrix(7,2)=(VF_Matrix(2,7)*A_1)/A_3; %%Aperature to panel 2
VF_Matrix(7,3)=(VF_Matrix(3,7)*A_1)/A_3; %%Aperature to panel 3
VF_Matrix(7,4)=(VF_Matrix(4,7)*A_1)/A_3; %%Aperature to panel 4
VF_Matrix(7,7)=0; %%Aperature to opening
VF_Matrix(7,5)=(1-sum(VF_Matrix(7,:)))/2; %%Aperature to floor
VF_Matrix(7,6)=VF_Matrix(7,5); %%Aperature to roof

%%%F-hat
KD=eye(N_el_domain);
F_hat_sol=inv(KD-VF_Matrix.*rho_sol')*VF_Matrix; %Solar F-Hat Matrix
F_hat_therm=inv(KD-VF_Matrix.*rho_therm')*VF_Matrix; %Thermal F-Hat Matrix

%%%Incoming Solar Flux
% Q_in=2.5e6; %(W), (kW e3),(MW e6)
Q_in=12.5e6; %(W), (kW e3),(MW e6)
Q_sol_flux(7,1)=0; %Initialize flux matrix
for i=1:4
    Q_sol_flux(i)=Q_in/A(1);
end

%%%Q_solar Radiation%%%
for i=1:7
    for j=1:7
        Q_sol_loss(i,j)=(Q_sol_flux(i)*A(i)*rho_sol(i))*F_hat_sol(i,j)*epsilon_sol(j)-(Q_sol_flux(j)*A(j)*rho_sol(j))*F_hat_sol(j,i)*epsilon_sol(i); %Reflected Solar energy
    end
end
Q_sol_loss_kW=sum(Q_sol_loss)/1000; %Solar radiation loss from each element
Q_in_kW=(Q_in*4)/1000; %Total energy incident to the system
Q_solar_loss_kW=Q_sol_loss_kW(7); %Solar energy lost from the aperature

%%%Form Ax=b, A-coeffecient matrix, x-temp vector, b-flux vector?
A_EB=zeros(6,6);
b_EB=zeros(6,1);

b_EB=(A.*Q_sol_flux')'+sum(Q_sol_loss)';
b_EB(7)=[];

for i=1:6
   b_EB(i)=b_EB(i)/(epsilon_therm(i)*A(i)*sigma)+epsilon_therm(7)*F_hat_therm(i,7)*T_infinity^4; 
end

for i=1:(N_el_domain-1)
    for j=1:(N_el_domain-1)
        A_EB(i,j)=-epsilon_therm(j)*F_hat_therm(i,j);
    end
    for j=1:(N_el_domain)
        A_EB(i,i)=A_EB(i,i)+epsilon_therm(j)*F_hat_therm(i,j);
    end
end
x=A_EB\b_EB;
T_max=x.^(1/4);
T=T_max;
% T_max(7)=T_infinity



% Avg_Temp=0;
% Cav_Area=0;
% for i=1:6
%    Cav_Area=Cav_Area+A(i);
% end
% for i=1:6
%    Avg_Temp=Avg_Temp+T(i)*(A(i)/Cav_Area); 
% end
% H=12; %Height of panel
% T_infinity=273.15+20;
% [h_bar_conv] = S_and_K_Function_v2(Avg_Temp, H, T_infinity);
h_bar_conv=0;

%%%Reset A and B vector
A_EB=zeros(6,6);
b_EB=zeros(6,1);
h_bar_conv=ones(6,1)*h_bar_conv;
UA_HTF=ones(6,1)*UA_HTF;
UA_HTF(5:6)=0;
clear T
factor=[0.9218; 0.9218; 0.9218; 0.9218; 0.6457; 0.6457];
T=T_max.*factor;
T(7)=T_infinity;
tol=1e-4;
error=tol*100;
iteration=1;

while tol < error
    A_EB=zeros(6,6);
    b_EB=zeros(6,1);
    for i=1:(N_el_domain-1)
        for j=1:(N_el_domain-1)
            A_EB(i,j)=-epsilon_therm(j)*F_hat_therm(i,j)*(T(i)^2+T(j)^2)*(T(i)+T(j));
        end
        for j=1:(N_el_domain)
            A_EB(i,i)=A_EB(i,i)+epsilon_therm(j)*F_hat_therm(i,j)*(T(i)^2+T(j)^2)*(T(i)+T(j));
        end
        A_EB(i,i)=A_EB(i,i)+(h_bar_conv(i)+UA_HTF(i))/(epsilon_therm(i)*sigma);
    end
    
    b_EB=(A.*Q_sol_flux')'+sum(Q_sol_loss)';
    b_EB(7)=[];
    for i=1:6
       b_EB(i)=b_EB(i)/(epsilon_therm(i)*A(i)*sigma)+epsilon_therm(7)*F_hat_therm(i,7)*T_infinity*(T(i)^2+T_infinity^2)*(T(i)+T_infinity)+(h_bar_conv(i)*T_infinity+UA_HTF(i)*T_avg_HTF)/(epsilon_therm(i)*sigma); 
    end
    
    T_star=A_EB\b_EB;
    T_star(7)=T_infinity;
    error=max(abs(T-T_star)./T);
    T=T_star;
%     Avg_Temp=0;
%     Cav_Area=0;
%     for i=1:6
%        Cav_Area=Cav_Area+A(i);
%     end
%     for i=1:6
%        Avg_Temp=Avg_Temp+T(i)*(A(i)/Cav_Area); 
%     end
%     H=12; %Height of panel
%     T_infinity=273.15+20;
%     [h_bar_conv] = S_and_K_Function_v2(Avg_Temp, H, T_infinity);
    h_bar_conv=0;
    h_bar_conv=ones(6,1)*h_bar_conv;
end
T;




for i=1:(N_el_domain-1)
   q_dot_conv_loss(i)=h_bar_conv(i)*A(i)*(T(i)-T_infinity); 
end

for i=1:4
   q_dot_HTF_gain(i)=UA_HTF(i)*A(i)*(T(i)-T_avg_HTF); 
end

for i=1:N_el_domain
    for j=1:N_el_domain
        q(i,j)=(epsilon_therm(i)*A(i)*sigma)*(epsilon_therm(j)*F_hat_therm(i,j)*(T(i)^4-T(j)^4));
    end
end
q_dot=sum(q,2);


q_conv_kW=q_dot_conv_loss/1000;
q_HTF_kW=q_dot_HTF_gain/1000;
q_dot_kW=(q_dot/1000)';

HTF_Gain_kW=sum(q_HTF_kW);
Conv_Loss_kW=sum(q_conv_kW);
Rad_Loss_kW=abs(q_dot_kW(7));


display('-----------------------------------------------')
display(['Energy In: ',num2str(Q_in_kW/1000), ' (MW)'])
display(['Solar Radation Loss: ',num2str(Q_solar_loss_kW/1000), ' (MW)'])
display(['Thermal Radation Loss: ',num2str(Rad_Loss_kW/1000), ' (MW)'])
display(['Convection Loss: ',num2str(Conv_Loss_kW/1000), ' (MW)'])
display(['HTF Gain: ',num2str(HTF_Gain_kW/1000), ' (MW)'])
display('-----------------------------------------------')
display(['Solar Radation Loss Percentage: ',num2str((Q_solar_loss_kW/Q_in_kW)*100), ' (%)'])
display(['Thermal Radation Loss Percentage: ',num2str((Rad_Loss_kW/Q_in_kW)*100), ' (%)'])
display(['Convection Loss Percentage: ',num2str((Conv_Loss_kW/Q_in_kW)*100), ' (%)'])
display(['HTF Gain Percentage: ',num2str((HTF_Gain_kW/Q_in_kW)*100), ' (%)'])
display(['Residual Percentage: ',num2str(100-((Rad_Loss_kW/Q_in_kW)*100+(Q_solar_loss_kW/Q_in_kW)*100+(Conv_Loss_kW/Q_in_kW)*100+(HTF_Gain_kW/Q_in_kW)*100)), ' (%)'])
display('-----------------------------------------------')
display(['Panel 1 Temperature: ' num2str(T(1)), ' (K)'])
display(['Panel 2 Temperature: ' num2str(T(2)), ' (K)'])
display(['Panel 3 Temperature: ' num2str(T(3)), ' (K)'])
display(['Panel 4 Temperature: ' num2str(T(4)), ' (K)'])
display(['Roof Temperature: ' num2str(T(5)), ' (K)'])
display(['Floor Temperature: ' num2str(T(6)), ' (K)'])
display(['Aperature: ' num2str(T(7)), ' (K)'])
display(['Average HTF Temp: ' num2str(T_avg_HTF), ' (K)'])
display('-----------------------------------------------')
