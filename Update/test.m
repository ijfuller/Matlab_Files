clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cavity Geometry %%%%%
Receiver_Height=12;  %Receiver opening height in meters
Receiver_Width=14;   %Reciever opening width in meters
Top_Lip=0;           %Height of top lip in meters
Bottom_Lip=0;        %Height of bottom lip in meters

%%%%% Cavity Mesh %%%%%
%%%Cover Mesh Options: cover_n==1 && cover_m==1; cover_n==2 && cover_m==1; cover_n==4 && cover_m==1; cover_n==2 && cover_m==2; cover_n==2 && cover_m==3; cover_n==4 && cover_m==3
panel_n=2; % Number of elements in x axis for panel
panel_m=1; % Number of elements in y axis for panel
cover_n=1; %Number of elements in x axis for top and bottom of cavity.
cover_m=1; %Number of elements in y axis for top and bottom of cavity.
Plot_Mesh=0; %Disable=0 Enable=1. This will plot the mesh of the cavity.

%%%%% Absorbtivity/Emissivity of Cavity Surfaces %%%%%
%%%Solar emissivity
e_act_sol=0.965; %Absorbtivity in short wave range for active surfaces
e_pass_sol=0.05;%Absorbtivity in short wave range for passive surfaces
%%%Thermal emissivity
e_act_therm=0.85; %Emissivity in long wave range for active surfaces
e_pass_therm=0.05; %Emissivity in long wave range for passive surfaces

%%%Incoming Solar%%%
Q_in=(100e6/4); %Incident solar radiation per panel
Dist_Pattern=1; %1-uniform distribution. 2-Distribution for elements in a single panel, where all panels have the same distribution. 3 - each panel and element distribution is defined
Q_nodal_perc=[0 0 0.1 0.2 0.2 0.2 0.2 0.1 0 0]; %For Dist_Pattern=2, need to define elmenets in 1 panel (Commas not semicolons, and percentage value in decimal form, ex 0.8 not 80%). For pattern 3, need to define every panel (x4) element decimal value.

%%%%% Heat Transfer Properties %%%%%
T_HTF_in=290+273.15;%Inlet HTF temperature in kelvin
T_HTF_out=575+273.15;%Outlet HTF temperature in kelvin
UA_Type=1; %1-Program solves for UA value. 2-User defined UA value. Overall heat transfer coeffeceint of piping to HTF
UA_HTF=4000; %UA value for user defined type.

%%%%% Air/Convection %%%%%
T_infinity=20+273.15;%Ambient air temperature in kelvin
h_Type=2; %Heat transfer coeffecient for convection. Type 1: User Defined, Type 2: Siebers & Kraabel, Type 3:Clausing 1987
h_bar_conv=0; %For user defined heat transfer coeffecient for convection

N_el_domain=(4*panel_n*panel_m)+(2*cover_n*cover_m)+1;
A_Panel=5.36*12;
A_Roof=(1+sqrt(2))*5.36^2;
A_Aperature=14*12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

epsilon_sol =[0.9650    0.9650    0.9650    0.9650    0.9650    0.9650    0.9650    0.9650    0.0500    0.0500    1.0000];


epsilon_therm =[0.8500    0.8500    0.8500    0.8500    0.8500    0.8500    0.8500    0.8500    0.0500    0.0500    1.0000];


Q_nodal_perc =[0.5000    0.5000    0.5000    0.5000    0.5000    0.5000    0.5000    0.5000];


F_hat_sol =[0.0153    0.0199    0.0381    0.0377    0.0582    0.0513    0.0635    0.0566    0.1387    0.1387    0.6575;
    0.0199    0.0261    0.0939    0.0454    0.0863    0.0633    0.0806    0.0635    0.1811    0.1811    0.5196;
    0.0381    0.0939    0.0252    0.0264    0.0450    0.0424    0.0633    0.0513    0.1788    0.1788    0.6104;
    0.0377    0.0454    0.0264    0.0277    0.0959    0.0450    0.0863    0.0582    0.1872    0.1872    0.5735;
    0.0582    0.0863    0.0450    0.0959    0.0277    0.0264    0.0454    0.0377    0.1872    0.1872    0.5735;
    0.0513    0.0633    0.0424    0.0450    0.0264    0.0252    0.0939    0.0381    0.1788    0.1788    0.6104;
    0.0635    0.0806    0.0633    0.0863    0.0454    0.0939    0.0261    0.0199    0.1811    0.1811    0.5196;
    0.0566    0.0635    0.0513    0.0582    0.0377    0.0381    0.0199    0.0153    0.1387    0.1387    0.6575;
    0.0643    0.0840    0.0829    0.0868    0.0868    0.0829    0.0840    0.0643    0.0162    0.1189    0.3795;
    0.0643    0.0840    0.0829    0.0868    0.0868    0.0829    0.0840    0.0643    0.1189    0.0162    0.3795;
    0.1258    0.0994    0.1168    0.1097    0.1097    0.1168    0.0994    0.1258    0.1565    0.1565    0.1125];


F_hat_therm =[0.0173    0.0227    0.0401    0.0401    0.0602    0.0533    0.0658    0.0581    0.1458    0.1458    0.6817;
    0.0227    0.0302    0.0965    0.0487    0.0890    0.0661    0.0838    0.0658    0.1912    0.1912    0.5536;
    0.0401    0.0965    0.0279    0.0290    0.0475    0.0448    0.0661    0.0533    0.1868    0.1868    0.6374;
    0.0401    0.0487    0.0290    0.0310    0.0983    0.0475    0.0890    0.0602    0.1960    0.1960    0.6032;
    0.0602    0.0890    0.0475    0.0983    0.0310    0.0290    0.0487    0.0401    0.1960    0.1960    0.6032;
    0.0533    0.0661    0.0448    0.0475    0.0290    0.0279    0.0965    0.0401    0.1868    0.1868    0.6374;
    0.0658    0.0838    0.0661    0.0890    0.0487    0.0965    0.0302    0.0227    0.1912    0.1912    0.5536;
    0.0581    0.0658    0.0533    0.0602    0.0401    0.0401    0.0227    0.0173    0.1458    0.1458    0.6817;
    0.0676    0.0887    0.0866    0.0909    0.0909    0.0866    0.0887    0.0676    0.0296    0.1322    0.4245;
    0.0676    0.0887    0.0866    0.0909    0.0909    0.0866    0.0887    0.0676    0.1322    0.0296    0.4245;
    0.1304    0.1059    0.1220    0.1154    0.1154    0.1220    0.1059    0.1304    0.1751    0.1751    0.1773];


rho_sol =[0.0350;
    0.0350;
    0.0350;
    0.0350;
    0.0350;
    0.0350;
    0.0350;
    0.0350;
    0.9500;
    0.9500;
         0];


rho_therm =[0.1500;
    0.1500;
    0.1500;
    0.1500;
    0.1500;
    0.1500;
    0.1500;
    0.1500;
    0.9500;
    0.9500;
         0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
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