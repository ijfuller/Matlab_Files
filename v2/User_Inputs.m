clear all
clc
close all

%%%User Inputs%%%
Receiver_Height=12;  %Receiver opening height
Receiver_Width=14;   %Reciever opening width
Top_Lip=2;
Bottom_Lip=2;

panel_n=1; % Number of elements in x axis for panel
panel_m=10; % Number of elements in y axis for panel

cover_n=1; %Number of elements in x axis for top and bottom of cavity.***Some issue, needs to be even***
cover_m=1; %Number of elements in y axis for top and bottom of cavity

%%%Solar emissivity
e_act_sol=0.95; %Emissivity in short wave range for active surfaces
e_pass_sol=0.47;%Emissivity in short wave range for passive surfaces
%%%Vector of solar emissivities
for i=1:4*panel_m
   epsilon_sol(i)=e_act_sol; 
end
for i=1+panel_m*4:2+panel_m*4
    epsilon_sol(i)=e_pass_sol;
end
epsilon_sol(3+panel_m*4)=1;

%%%Thermal emissivity
e_act_therm=0.95; %Emissivity in long wave range for active surfaces
e_pass_therm=0.47; %Emissivity in long wave range for passive surfaces
%%%Vector of thermal emissivities
for i=1:4*panel_m
   epsilon_therm(i)=e_act_therm; 
end
for i=1+panel_m*4:2+panel_m*4
    epsilon_therm(i)=e_pass_therm;
end
epsilon_therm(3+panel_m*4)=1;


Q_in=(250e6/4);%50e6/4; %Incident solar radiation per panel
T_HTF_in=290+273.15;%In kelvin
T_HTF_out=575+273.15;%In kelvin
T_infinity=20+273.15;
UA_HTF=2000;
h_Type=2; %Type 1: User Defined, Type 2: Siebers & Kraabel, Type 3:Clausing 1987
h_bar_conv=8.6098; %For user deefined heat transfer coeffecient


Dist_Pattern=2; %1 is uniform, 2 is panel uniform (needs 1 set of percentages in decimal format), 3 differnt panels
Q_nodal_perc=[0 0 0 0 0 0 0 0 .5 .5]; %(Commas not semicolons)
%Q_nodal_perc=[0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %(Commas not semicolons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Dist_Pattern==1
    Q_nodal_perc=ones(1,panel_m*4)*(1/panel_m);
elseif Dist_Pattern==2
    [i,j]=size(Q_nodal_perc);
    if j > panel_m || j < panel_m
        display("Error! Q_nodal_perc dimensions do not match number of nodes.")
    end
    temp=Q_nodal_perc;
    for i=2:4
        temp=[temp Q_nodal_perc]; 
    end
    Q_nodal_perc=temp;
else 
    [i,j]=size(Q_nodal_perc);
    if j > panel_m*4 || j < panel_m*4
        display("Error! Q_nodal_perc dimensions do not match number of nodes.")
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cal_View_Factors=0; %Disable=0,Enable=1;
N_rays=5e6; %Number of rays for Monte Carlo simulation
Save_View_Factor_File=1;%Disable=0 Enable=1. Will export view factor data.
Save_View_Factor_Filename=['VF_Matrix_',num2str(panel_m),'_Nodes_Lip_2m.xlsx']; %Will only save the file if Save_View_Factor_File=1.

%%% OR Load File is Cal_Veiw Factors=0 %%%
Load_View_Factor_Matrix_File=['VF_Matrix_',num2str(panel_m),'_Nodes_Lip_2m.xlsx']; %Will only load the file if Cal_View_Factors=0.

%%%Plotting%%%
Plot_Mesh=0; %Disable=0 Enable=1. This will plot the mesh of the cavity. (Seperatly)


%%%%Going to move this all to a subfunction routine%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Code Portion%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Call cavity mesh function
[El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, N_el_domain, A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip] = Cavity_Mesh_Function(Receiver_Height,Receiver_Width,panel_n,panel_m,cover_n,cover_m,Plot_Mesh,Top_Lip,Bottom_Lip);

if Cal_View_Factors==1
    %%%Call F3D_30 View Factor (Panel to Panel view factors)
    %F3D_30_Function(panel_m,panel_n,Receiver_Width,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel); %Not sure why this is commented out
    [F_Panel_1_to_Panel, F_Panel_2_to_Panel] = F3D_30_Function_v3(panel_m,panel_n,Receiver_Width,Receiver_Height,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel,Top_Lip,Bottom_Lip);

    %%%Call Monte Carlo Function - Panel 1 & 2 to floor 
    if cover_n==1 && cover_m==1
        [F_Panel_1_Floor  F_Panel_2_Floor F_Aper_Floor,F_Aper_Roof, F_Bottom_Lip_Floor, F_Top_Lip_Floor,  F_Bottom_Lip_Roof, F_Top_Lip_Roof] = Monte_Carlo_Panel_Floor_Pent_Function_v4(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, Receiver_Width, Receiver_Height,Top_Lip, Bottom_Lip);
    else
        %[F_Panel_1_Floor] = Monte_Carlo_Panel_Floor_Function(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover);
    end

    %%%Call Monte_Carlo Function Floor to Roof
    if cover_n==1 && cover_m==1
        [F_Roof_Floor] = Monte_Carlo_roof_to_floor_Function(Receiver_Height,Receiver_Width,N_rays,Global_Nodal_Coor_Cover);
    end

    %%%Build View Factor Matrix
    [VF_Matrix] = View_Factor_Matrix_Function_v3(N_el_domain,F_Panel_1_to_Panel,F_Panel_2_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename,panel_m);
else
     VF_Matrix=xlsread(Load_View_Factor_Matrix_File);
end

%%%F Hat Matrix%%%
[F_hat_sol,F_hat_therm,epsilon_sol,epsilon_therm,rho_sol,rho_therm] = F_Hat_Function_v3(N_el_domain,VF_Matrix, epsilon_sol,epsilon_therm,Top_Lip,Bottom_Lip,e_pass_sol,e_pass_therm,panel_m);

%%%Energy Balance%%%
[T_max,T,Q_in_kW,Q_solar_loss_kW,Rad_Loss_kW,Conv_Loss_kW,HTF_Gain_kW,T_Nodal] = Energy_Balance_Function_v5(N_el_domain,epsilon_sol,epsilon_therm,rho_sol,rho_therm,F_hat_sol,F_hat_therm,T_HTF_in,T_HTF_out,T_infinity,A_Aperature,A_Panel,A_Roof,Receiver_Height,A_Top_Lip,A_Bottom_Lip,UA_HTF,h_Type,h_bar_conv,Q_in,panel_m,Top_Lip,Bottom_Lip,Q_nodal_perc);

% %A_Total=(A_Roof*2+A_Top_Lip*2);
% %Av_Temp_Passive=(T(5)*A_Roof/A_Total)+(T(5)*A_Roof/A_Total)+(T(7)*A_Top_Lip/A_Total)+(T(7)*A_Top_Lip/A_Total)
