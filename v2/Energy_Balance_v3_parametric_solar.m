clear all
clc

e_act_sol=[0:0.01:1];
e_pass_sol=[0:0.2:1];

for ii=1:length(e_act_sol)
    for jj=1:length(e_pass_sol)
        %%%Solar emissivity
%         e_act_sol=0.5;
%         e_pass_sol=1;
        epsilon_sol=[e_act_sol(ii);e_act_sol(ii);e_act_sol(ii);e_act_sol(ii);e_pass_sol(jj);e_pass_sol(jj);1];
        rho_sol=1-epsilon_sol;
        %%%Thermal emissivity
        e_act_therm=1;
        e_pass_therm=1;
        epsilon_therm=[e_act_therm;e_act_therm;e_act_therm;e_act_therm;e_pass_therm;e_pass_therm;1];
        rho_therm=1-epsilon_therm;

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
        F_hat_sol=inv(KD-VF_Matrix.*rho_sol')*VF_Matrix;
        F_hat_therm=inv(KD-VF_Matrix.*rho_therm')*VF_Matrix;

        %%%Incoming Solar Flux
        Q_in=100000; %(W)
        Q_sol_flux(7,1)=0;
        for i=1:4
            Q_sol_flux(i)=Q_in/A(1);
        end

        %%%Q_solar Radiation%%%
        for i=1:7
            for j=1:7
                Q_sol_loss(i,j)=(Q_sol_flux(i)*A(i)*rho_sol(i))*F_hat_sol(i,j)*epsilon_sol(j)-(Q_sol_flux(j)*A(j)*rho_sol(j))*F_hat_sol(j,i)*epsilon_sol(i);
            end
        end
        Q_sol_loss_kW=sum(Q_sol_loss)/1000;

        Q_in_kW=(Q_in*4)/1000;
        Q_solar_loss_kW=Q_sol_loss_kW(7);
        Q_sol_loss_perc(ii,jj)=(Q_solar_loss_kW/Q_in_kW)*100;
    end
end
Q_sol_loss_perc

figure(1)
hold on
grid on
for jj=1:length(e_pass_sol)
    plot(e_act_sol,Q_sol_loss_perc(:,jj),'Linewidth',2)
end
title('Q Solar Loss %')
xlabel('Solar Epsilon Active') 
ylabel('Loss %') 
legend({'e_p_a_s_s_i_v_e = 0','e_p_a_s_s_i_v_e = 0.2','e_p_a_s_s_i_v_e = 0.4','e_p_a_s_s_i_v_e = 0.6','e_p_a_s_s_i_v_e = 0.8','e_p_a_s_s_i_v_e = 1'},'Location','southwest')