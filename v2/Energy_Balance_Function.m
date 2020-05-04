function [] = Energy_Balance_Function(N_el_domain,epsilon_sol,F_Hat,T_HTF_in,T_HTF_out,T_infinity,A_Aperature,A_Panel,A_Roof,Receiver_Height)

    F_hat=F_Hat;
    epsilon=epsilon_sol;
    rho=1-epsilon;
    UA_HTF=150;
    T_avg_HTF=(T_HTF_in+T_HTF_out)/2;
    sigma=5.67*10^-8;


    A_1=A_Panel;%Panel_Area
    A_2=A_Roof; %Roof Area
    A_3=A_Aperature; %Aperature Area

    A(1)=A_1;
    A(2)=A_1;
    A(3)=A_1;
    A(4)=A_1;
    A(5)=A_2;
    A(6)=A_2;
    A(7)=A_3;

    Q_solar=1000*1000*10;
    %%%Form Ax=b, A-coeffecient matrix, x-temp vector, b-
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

    %Getting starting point for h_bar function
    Avg_Temp=0;
    Cav_Area=0;
    for i=1:6
       Avg_Temp=Avg_Temp+T_max(i)*A(i);
       Cav_Area=Cav_Area+A(i);
    end
    H=Receiver_Height; %Height of panel
    Avg_Temp=Avg_Temp/Cav_Area;
    [h_bar_conv] = S_and_K_Function(Avg_Temp, H, T_infinity);

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
        error=max(abs(T-T_star)./T);
        T=T_star;
        for i=1:6
            Avg_Temp=Avg_Temp+T(i)*A(i);
        end
        Avg_Temp=Avg_Temp/Cav_Area;
        [h_bar_conv] = S_and_K_Function(Avg_Temp, H, T_infinity);
        h_bar_conv=ones(6,1)*h_bar_conv;
        iteration=iteration+1;
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
            q(i,j)=(epsilon(i)*A(i)*sigma)*(epsilon(j)*F_hat(i,j)*(T(i)^4-T(j)^4));
        end
    end
    q_dot=sum(q,2);

    q_conv_kW=q_dot_conv_loss/1000;
    q_HTF_kW=q_dot_HTF_gain/1000;
    q_dot_kW=(q_dot/1000)';

    Q_in_kW=(Q_solar*4/1000);
    HTF_Gain_kW=sum(q_HTF_kW);
    Conv_Loss_kW=sum(q_conv_kW);
    Rad_Loss_kW=abs(q_dot_kW(7));
    Res_Perc=((Q_in_kW-(HTF_Gain_kW+Conv_Loss_kW+Rad_Loss_kW))/Q_in_kW)*100;


    display('----------------------------------------')
    display(['Energy In: ',num2str(Q_in_kW), ' (kW)'])
    display(['Radation Loss: ',num2str(Rad_Loss_kW), ' (kW)'])
    display(['Convection Loss: ',num2str(Conv_Loss_kW), ' (kW)'])
    display(['HTF Gain: ',num2str(HTF_Gain_kW), ' (kW)'])
    display('----------------------------------------')
    display(['Radation Loss Percentage: ',num2str((Rad_Loss_kW/Q_in_kW)*100), ' (%)'])
    display(['Convection Loss Percentage: ',num2str((Conv_Loss_kW/Q_in_kW)*100), ' (%)'])
    display(['HTF Gain Percentage: ',num2str((HTF_Gain_kW/Q_in_kW)*100), ' (%)'])
    display(['Residual Percentage: ',num2str(100-((Rad_Loss_kW/Q_in_kW)*100+(Conv_Loss_kW/Q_in_kW)*100+(HTF_Gain_kW/Q_in_kW)*100)), ' (%)'])
    display('----------------------------------------')
    display(['Panel 1 Temperature: ' num2str(T(1)), ' (K)'])
    display(['Panel 2 Temperature: ' num2str(T(2)), ' (K)'])
    display(['Panel 3 Temperature: ' num2str(T(3)), ' (K)'])
    display(['Panel 4 Temperature: ' num2str(T(4)), ' (K)'])
    display(['Roof Temperature: ' num2str(T(5)), ' (K)'])
    display(['Floor Temperature: ' num2str(T(6)), ' (K)'])
    display(['Aperature: ' num2str(T(7)), ' (K)'])
    display('----------------------------------------')
end