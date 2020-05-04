
clear all
clc

% e_active=0.9999
active=[0.01:0.01:0.99];
for par=1:length(active)
e_active=active(par)
e_passive=0.5
epsilon=[e_active;e_active;e_active;e_active;e_passive;e_passive;1];
rho=1-epsilon;
%UA_HTF=15;
UA_HTF=100;
T_infinity=273.15+20;
T_avg_HTF=400+273.15;


N_el_domain=7;
VF_Matrix=zeros(N_el_domain,N_el_domain); %Initialize fiew factor full matrix
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

%E=(Temp.^4).*(5.67*10^-8);

KD=eye(N_el_domain);
VF_Matrix;
F_hat=inv(KD-VF_Matrix.*rho')*VF_Matrix;

% %%%Convert to Standard Form%%%
% N_surf=7;
% F=reshape(VF_Matrix',[N_surf^2,1]);
% A=zeros(N_surf^2,N_surf^2);
% 
% int=1;
% for i=1:N_surf
%     for j=1:N_surf
%         for k=1:N_surf
%             j+(i-1)*N_surf;
%             Temp(int)=F(k+(i-1)*N_surf);
%             int=1+int;
%         end
%     end
% end
% Temp=reshape(Temp,[N_surf,N_surf^2])';
% Temp=Temp.*rho';
% 
% k=1;
% for i=1:N_surf^2
%         Temp(i,k)=Temp(i,k)-1;
%     if mod(i,N_surf)==0
%         k=k+1;
%     end
% end
% 
% for i=1:N_surf^2
%     Set(i,:)=[1:N_surf:N_surf^2];
% end
% for i=2:N_surf
%     Rows=[i:N_surf:N_surf^2];
%     for j=1:N_surf
%         Set(Rows(j),:)=Set(Rows(j),:)+i-1;
%     end
% end
% 
% A=zeros(N_surf^2,N_surf^2);
% for i=1:N_surf^2
%     for j=1:N_surf
%         A(i,Set(i,j))=Temp(i,j); 
%     end
% end
% A;
%A\F


    Q_solar=1000*1000*10;
    %%%Form Ax=b, A-coeffecient matrix, x-temp vector, b-flux vector?
    A_EB=zeros(6,6);
    b_EB=zeros(6,1);

    for i=1:4
        b_EB(i)=Q_solar;
    end
    for i=1:6
       b_EB(i)=b_EB(i)/(epsilon(i)*A(i)*sigma)+epsilon(7)*F_hat(i,7)*T_infinity^4; 
    end

    for i=1:(N_el_domain-1)
        for j=1:(N_el_domain-1)
            A_EB(i,j)=-epsilon(j)*F_hat(i,j);
        end
        for j=1:(N_el_domain)
            A_EB(i,i)=A_EB(i,i)+epsilon(j)*F_hat(i,j);
        end
    end
    x=A_EB\b_EB;
    T_max=x.^(1/4);

    T=T_max;
    T(7)=T_infinity

    %Getting starting point for h_bar function
    Avg_Temp=0;
    Cav_Area=0;
    for i=1:6
       Avg_Temp=Avg_Temp+T_max(i)*A(i);
       Cav_Area=Cav_Area+A(i);
    end
    H=12; %Height of panel
    Avg_Temp=Avg_Temp/Cav_Area
    T_infinity=273.15;
    [h_bar_conv] = S_and_K_Function(Avg_Temp, H, T_infinity)

    %%%Reset A and B vector
    A_EB=zeros(6,6);
    b_EB=zeros(6,1);
    h_bar_conv=ones(6,1)*h_bar_conv;
    UA_HTF=ones(6,1)*UA_HTF;
    UA_HTF(5:6)=0;
    clear T
    factor=[0.9218; 0.9218; 0.9218; 0.9218; 0.6457; 0.6457];
    T=T_max.*factor
    T(7)=T_infinity;
    tol=1e-4;
    error=tol*100;
    iteration=1

    while tol < error
        A_EB=zeros(6,6);
        b_EB=zeros(6,1);
        for i=1:(N_el_domain-1)
            for j=1:(N_el_domain-1)
                A_EB(i,j)=-epsilon(j)*F_hat(i,j)*(T(i)^2+T(j)^2)*(T(i)+T(j));
            end
            for j=1:(N_el_domain)
                A_EB(i,i)=A_EB(i,i)+epsilon(j)*F_hat(i,j)*(T(i)^2+T(j)^2)*(T(i)+T(j));
            end
            A_EB(i,i)=A_EB(i,i)+(h_bar_conv(i)+UA_HTF(i))/(epsilon(i)*sigma);
        end

        for i=1:4
            b_EB(i)=Q_solar;
        end
        for i=1:6
           b_EB(i)=b_EB(i)/(epsilon(i)*A(i)*sigma)+epsilon(7)*F_hat(i,7)*T_infinity*(T(i)^2+T_infinity^2)*(T(i)+T_infinity)+(h_bar_conv(i)*T_infinity+UA_HTF(i)*T_avg_HTF)/(epsilon(i)*sigma); 
        end

        T_star=A_EB\b_EB;
        T_star(7)=T_infinity;
        error=max(abs(T-T_star)./T)
        T=T_star
        for i=1:6
            Avg_Temp=Avg_Temp+T(i)*A(i);
        end
        Avg_Temp=Avg_Temp/Cav_Area;
        [h_bar_conv] = S_and_K_Function(Avg_Temp, H, T_infinity)
        h_bar_conv=ones(6,1)*h_bar_conv;
        iteration=iteration+1;
    end
    T




    for i=1:(N_el_domain-1)
       q_dot_conv_loss(i)=h_bar_conv(i)*A(i)*(T(i)-T_infinity); 
    end

    for i=1:4
       q_dot_HTF_gain(i)=UA_HTF(i)*A(i)*(T(i)-T_avg_HTF); 
    end

    for i=1:N_el_domain
        for j=1:N_el_domain
            q(i,j)=(epsilon(i)*A(i)*sigma)*(epsilon(j)*F_hat(i,j)*(T(i)^4-T(j)^4));
        end
    end
    q_dot=sum(q,2);

    q_conv_kW=q_dot_conv_loss/1000;
    q_HTF_kW=q_dot_HTF_gain/1000;
    q_dot_kW=(q_dot/1000)';

    Q_in_kW(par)=(Q_solar*4/1000);
    HTF_Gain_kW(par)=sum(q_HTF_kW);
    Conv_Loss_kW(par)=sum(q_conv_kW);
    Rad_Loss_kW(par)=abs(q_dot_kW(7));
    
    HTF_Gain_Percent(par)=(HTF_Gain_kW/Q_in_kW)*100;
    Conv_Loss_Percent(par)=(Conv_Loss_kW/Q_in_kW)*100;
    Rad_Loss_Percent(par)=(Rad_Loss_kW/Q_in_kW)*100;
    
    Res_Perc=((Q_in_kW-(HTF_Gain_kW+Conv_Loss_kW+Rad_Loss_kW))/Q_in_kW)*100;


%     display('----------------------------------------')
%     display(['Energy In: ',num2str(Q_in_kW), ' (kW)'])
%     display(['Radation Loss: ',num2str(Rad_Loss_kW), ' (kW)'])
%     display(['Convection Loss: ',num2str(Conv_Loss_kW), ' (kW)'])
%     display(['HTF Gain: ',num2str(HTF_Gain_kW), ' (kW)'])
%     display('----------------------------------------')
%     display(['Radation Loss Percentage: ',num2str((Rad_Loss_kW/Q_in_kW)*100), ' (%)'])
%     display(['Convection Loss Percentage: ',num2str((Conv_Loss_kW/Q_in_kW)*100), ' (%)'])
%     display(['HTF Gain Percentage: ',num2str((HTF_Gain_kW/Q_in_kW)*100), ' (%)'])
%     display(['Residual Percentage: ',num2str(100-((Rad_Loss_kW/Q_in_kW)*100+(Conv_Loss_kW/Q_in_kW)*100+(HTF_Gain_kW/Q_in_kW)*100)), ' (%)'])
%     display('----------------------------------------')
%     display(['Panel 1 Temperature: ' num2str(T(1)), ' (K)'])
%     display(['Panel 2 Temperature: ' num2str(T(2)), ' (K)'])
%     display(['Panel 3 Temperature: ' num2str(T(3)), ' (K)'])
%     display(['Panel 4 Temperature: ' num2str(T(4)), ' (K)'])
%     display(['Roof Temperature: ' num2str(T(5)), ' (K)'])
%     display(['Floor Temperature: ' num2str(T(6)), ' (K)'])
%     display(['Aperature: ' num2str(T(7)), ' (K)'])
%     display('----------------------------------------')
    
end

figure (1)
hold on
grid on
plot(active,HTF_Gain_Percent,'LineWidth',2)
plot(active,Conv_Loss_Percent,'LineWidth',2)
plot(active,Rad_Loss_Percent,'LineWidth',2)
legend('HTF Gain (%)','Convection Loss (%)', 'Radiation Loss (%)')