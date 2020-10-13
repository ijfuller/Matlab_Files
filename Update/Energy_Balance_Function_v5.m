function [T_max,T,Q_in_kW,Q_solar_loss_kW,Rad_Loss_kW,Conv_Loss_kW,HTF_Gain_kW,T_Nodal,h_bar_conv] = Energy_Balance_Function_v5(N_el_domain,epsilon_sol,epsilon_therm,rho_sol,rho_therm,F_hat_sol,F_hat_therm,T_HTF_in,T_HTF_out,T_infinity,A_Aperature,A_Panel,A_Roof,Receiver_Height,A_Top_Lip,A_Bottom_Lip,UA_HTF,h_Type,h_bar_conv,Q_in,panel_n,panel_m,cover_n,cover_m,Top_Lip,Bottom_Lip,Q_nodal_perc)
    T_max=0;
    Q_in_kW=0;
    T=0;
    Q_solar_loss_kW=0;
    Rad_Loss_kW=0;
    Conv_Loss_kW=0;
    HTF_Gain_kW=0;
    
    sigma=5.67*10^-8;
    T_avg_HTF=(T_HTF_out+T_HTF_in)/2;
    T_avg_HTF_bounds=[T_HTF_in:((T_HTF_out-T_HTF_in)/panel_m):T_HTF_out];
    for i=1:panel_m
        T_avg_HTF_nodal(i)=(T_avg_HTF_bounds(i)+T_avg_HTF_bounds(i+1))/2;
    end
    T_avg_HTF_nodal=T_avg_HTF_nodal.*ones(panel_n,1);
    T_avg_HTF_nodal=reshape(T_avg_HTF_nodal,[panel_n*panel_m,1]);
    T_avg_HTF_nodal=T_avg_HTF_nodal.*ones(1,4);
    T_avg_HTF_nodal=reshape(T_avg_HTF_nodal,[4*panel_n*panel_m,1]);
    T_avg_HTF=zeros(N_el_domain,1);
    for ii=1:(4*panel_n*panel_m)
      T_avg_HTF(ii)=T_avg_HTF_nodal(ii);
    end
    
    for i=1:(4*panel_n*panel_m)
        A(i)=A_Panel/(panel_n*panel_m);
    end
    for i=(1+(panel_n*panel_m*4)):((2*cover_n*cover_m)+(panel_n*panel_m*4))
        A(i)=A_Roof;
    end
    A(1+(4*panel_n*panel_m)+(2*cover_n*cover_m))=A_Aperature;
    
    if Top_Lip > 0
        A(1+(4*panel_n*panel_m)+(2*cover_n*cover_m))=A_Top_Lip;
        A(2+(4*panel_n*panel_m)+(2*cover_n*cover_m))=A_Aperature;
    end
    if Bottom_Lip > 0
        A(2+(4*panel_n*panel_m)+(2*cover_n*cover_m))=A_Bottom_Lip;
        A(3+(4*panel_n*panel_m)+(2*cover_n*cover_m))=A_Aperature;
    end
    
    %%%Incoming Solar Flux
    Q_sol_flux(N_el_domain,1)=0; %Initialize flux matrix
    Q_nodal=Q_nodal_perc*Q_in;
    for i=1:4*panel_n*panel_m
        Q_sol_flux(i)=Q_nodal(i)/(A_Panel/(panel_n*panel_m));
    end
    
    
    %%%Q_solar Radiation%%%
    for i=1:N_el_domain
        for j=1:N_el_domain
            Q_sol_loss(i,j)=(Q_sol_flux(i)*A(i)*rho_sol(i))*F_hat_sol(i,j)*epsilon_sol(j)-(Q_sol_flux(j)*A(j)*rho_sol(j))*F_hat_sol(j,i)*epsilon_sol(i); %Reflected Solar energy
        end
    end
    Q_sol_loss_kW=sum(Q_sol_loss)/1000; %Solar radiation loss from each element
    Q_in_kW=(Q_in*4)/1000; %Total energy incident to the system
    Q_solar_loss_kW=Q_sol_loss_kW(N_el_domain); %Solar energy lost from the aperature
    
    %%%Form Ax=b, A-coeffecient matrix, x-temp vector, b-flux vector?
    A_EB=zeros((N_el_domain-1),(N_el_domain-1));
    b_EB=zeros((N_el_domain-1),1);

    b_EB=(A.*Q_sol_flux')'+sum(Q_sol_loss)';
    b_EB(N_el_domain)=[];
    
    for i=1:(N_el_domain-1)
        b_EB(i)=b_EB(i)/(epsilon_therm(i)*A(i)*sigma)+epsilon_therm(N_el_domain)*F_hat_therm(i,N_el_domain)*T_infinity^4; 
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
    T_h=T;
    if h_Type==2
        [h_bar_conv] = H_Correlations_Function(Receiver_Height,T_h,A,T_infinity,h_Type);
    end
    %%%Reset A and B vector
    A_EB=zeros((N_el_domain-1),(N_el_domain-1));
    b_EB=zeros((N_el_domain-1),1);
    h_bar_conv=ones((N_el_domain-1),1)*h_bar_conv;
    UA_HTF=ones((N_el_domain-1),1)*UA_HTF;
    UA_HTF((1+(panel_n*panel_m*4)):N_el_domain)=0;
    clear T
    T=T_max;
    T(N_el_domain)=T_infinity;
    tol=1e-4;
    error=tol*100;
    iteration=1;
    while tol < error
        A_EB=zeros((N_el_domain-1),(N_el_domain-1));
        b_EB=zeros((N_el_domain-1),1);
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
        b_EB(N_el_domain)=[];

        for i=1:(N_el_domain-1)
           b_EB(i)=b_EB(i)/(epsilon_therm(i)*A(i)*sigma)+epsilon_therm(N_el_domain)*F_hat_therm(i,N_el_domain)*T_infinity*(T(i)^2+T_infinity^2)*(T(i)+T_infinity)+(h_bar_conv(i)*T_infinity+UA_HTF(i)*T_avg_HTF(i))/(epsilon_therm(i)*sigma); 
        end

        T_star=A_EB\b_EB;
        T_star(N_el_domain)=T_infinity;
        error=max(abs(T-T_star)./T);
        T=T_star;
        if h_Type==2
            T_h=T;
            T_h(length(T_h))=[];
            [h_bar_conv] = H_Correlations_Function(Receiver_Height,T,A,T_infinity,h_Type);
        end
        if length(h_bar_conv)==1
            h_bar_conv=ones((N_el_domain-1),1)*h_bar_conv;
        else
            
        end
        iteration=iteration+1;
    end

    for i=1:(N_el_domain-1)
       q_dot_conv_loss(i)=h_bar_conv(i)*A(i)*(T(i)-T_infinity); 
    end

    for i=1:(4*panel_m*panel_n)
       q_dot_HTF_gain(i)=UA_HTF(i)*A(i)*(T(i)-T_avg_HTF(i)); 
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
    Rad_Loss_kW=abs(q_dot_kW(N_el_domain));
     
    %%%Recombine%%%
    T_Nodal=T;
    if panel_m*panel_n > 1
        T=zeros(7,1);
        for i=1:4
            for j=1:panel_m*panel_n
               T(i)=T(i)+T_Nodal(j+((i-1)*panel_m*panel_n))/(panel_m*panel_n);
            end
        end
        T(5)=T_Nodal(1+panel_m*4);
        T(6)=T_Nodal(2+panel_m*4);
        T(7)=T_infinity;
        if Top_Lip > 0
            T(7)=T_Nodal(3+panel_m*panel_n*4);
            T(8)=T_infinity;
        end
        if Bottom_Lip > 0
            T(8)=T_Nodal(4+panel_m*panel_n*4);
            T(9)=T_infinity;
        end
    end

    
    if h_Type==1
        display('User Defined h value')
        display(['Heat Transfer Coeffecient: ',num2str(h_bar_conv(1)), ' (W/m^2)'])
    else
        display('Siebers & Kraabel Correlation')
        display(['Heat Transfer Coeffecient: ',num2str(h_bar_conv(1)), ' (W/m^2)'])
    end
%     display('-----------------------------------------------')
%     display(['Energy In: ',num2str(Q_in_kW), ' (kW)'])
%     display(['Solar Radation Loss: ',num2str(Q_solar_loss_kW), ' (kW)'])
%     display(['Thermal Radation Loss: ',num2str(Rad_Loss_kW), ' (kW)'])
%     display(['Convection Loss: ',num2str(Conv_Loss_kW), ' (kW)'])
%     display(['HTF Gain: ',num2str(HTF_Gain_kW), ' (kW)'])
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
    if Top_Lip > 0 && Bottom_Lip==0
        display(['Top Lip: ' num2str(T(7)), ' (K)'])
        display(['Aperature: ' num2str(T(8)), ' (K)'])
    elseif Top_Lip > 0 && Bottom_Lip > 0
        display(['Top Lip: ' num2str(T(7)), ' (K)'])
        display(['Bottom Lip: ' num2str(T(8)), ' (K)'])
        display(['Aperature: ' num2str(T(9)), ' (K)'])
    else
        display(['Aperature: ' num2str(T(7)), ' (K)'])
    end
    display('-----------------------------------------------')
end