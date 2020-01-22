clear all
clc
close all
% tic
%%%User Inputs%%%
Receiver_Height=12;  %Receiver opening height
Receiver_Width=14;   %Reciever opening width

panel_n=1; % Number of elements in x axis for panel
panel_m=1; % Number of elements in y axis for panel

cover_n=1; %Number of elements in x axis for top and bottom of cavity. %Some issue, needs to be even
cover_m=1; %Number of elements in y axis for top and bottom of cavity


    N_rays=1e7; %Number of rays for random

    %%%%%%%%%%%Code Portion%%%%%%%%%%%%%%%%%%%%%%%

    %%%Call cavity mesh function
    [El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, N_el_domain]=Cavity_Mesh_Function(Receiver_Height,Receiver_Width,panel_n,panel_m,cover_n,cover_m);

    %%%Call F3D_30 View Factor (Panel to Panel view factors)
    %[F_Panel_to_Panel] = F3D_30_Function(panel_m,panel_n,Receiver_Width,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel);
    [F_Panel_1_to_Panel, F_Panel_2_to_Panel] = F3D_30_Function_v3(panel_m,panel_n,Receiver_Width,Receiver_Height,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel);

    %%%Call Monte Carlo Function - Panel 1 & 2 to floor 
    if cover_n==1 && cover_m==1
        [F_Panel_1_Floor  F_Panel_2_Floor F_Aper_Floor] = Monte_Carlo_Panel_Floor_Pent_Function_v3(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, Receiver_Width, Receiver_Height);
    else
        [F_Panel_1_Floor] = Monte_Carlo_Panel_Floor_Function(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover);
    end

    %%%Call Monte_Carlo Function Floor to Roof
%     (F_Panel_to_Panel(:,:,1)+F_Panel_to_Panel(:,:,2)+F_Panel_to_Panel(:,:,3)+F_Panel_to_Panel(:,:,4)+F_Panel_1_Floor*2)
