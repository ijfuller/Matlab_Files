clear all
clc
close all

%%%%% User Inputs for Simulation %%%%%

%%%%% Cavity Geometry %%%%%
Receiver_Height=12;  %Receiver opening height in meters
Receiver_Width=14;   %Reciever opening width in meters
Top_Lip=0;           %Height of top lip in meters
Bottom_Lip=0;        %Height of bottom lip in meters

%%%%% Cavity Mesh %%%%%
%%%Cover Mesh Options: cover_n==1 && cover_m==1; cover_n==2 && cover_m==1; cover_n==4 && cover_m==1; cover_n==2 && cover_m==2; cover_n==2 && cover_m==3; cover_n==4 && cover_m==3
panel_n=1; % Number of elements in x axis for panel
panel_m=1; % Number of elements in y axis for panel
cover_n=1; %%%TEMPORARY, DO NOT CHANGE    %Number of elements in x axis for top and bottom of cavity.
cover_m=1; %%%TEMPORARY, DO NOT CHANGE    %Number of elements in y axis for top and bottom of cavity.
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

%%%%% Monte Carlo Simulation %%%%%
Cal_View_Factors=1; %Disable=0,Enable=1. Will calulate view factors if enabled.;
N_rays=1e7; %Number of rays for Monte Carlo simulation
Save_View_Factor_File=1;%Disable=0 Enable=1. Will export view factor data to excel.
Save_View_Factor_Filename=['VF_Matrix_Panel_n_',num2str(panel_n),'_Panel_m_',num2str(panel_m),'_Cover_n_',num2str(cover_n),'_Cover_m_',num2str(cover_m),'_Nodes.xlsx']; %View factor file name. Will only save the viewfactors if Save_View_Factor_File=1.
Load_View_Factor_Matrix_File=['VF_Matrix_Panel_n_',num2str(panel_n),'_Panel_m_',num2str(panel_m),'_Cover_n_',num2str(cover_n),'_Cover_m_',num2str(cover_m),'_Nodes.xlsx']; %Load view factors file. Will only load the file if Cal_View_Factors=0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function will send user defined variables to Function_Set to begin the simulation %%%
Function_Set(Receiver_Height,Receiver_Width,Top_Lip,Bottom_Lip,panel_m,panel_n,cover_n,cover_m,e_act_sol,e_pass_sol,e_act_therm,e_pass_therm,Q_in,T_HTF_in,T_HTF_out,T_infinity,UA_Type,UA_HTF,h_Type,h_bar_conv,Dist_Pattern,Q_nodal_perc,Cal_View_Factors,N_rays,Save_View_Factor_File,Save_View_Factor_Filename,Load_View_Factor_Matrix_File,Plot_Mesh)
