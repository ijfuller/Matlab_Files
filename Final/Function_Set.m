function [] = Function_Set(Receiver_Height,Receiver_Width,Top_Lip,Bottom_Lip,panel_m,panel_n,cover_n,cover_m,e_act_sol,e_pass_sol,e_act_therm,e_pass_therm,Q_in,T_HTF_in,T_HTF_out,T_infinity,UA_Type,UA_HTF,h_Type,h_bar_conv,Dist_Pattern,Q_nodal_perc,Cal_View_Factors,N_rays,Save_View_Factor_File,Save_View_Factor_Filename,Load_View_Factor_Matrix_File,Plot_Mesh);
%%%%% This function contains all the functions for the cavity simulation %%%%%

%%%%% Mesh Function: Meshes Panels and Roof/Floor %%%%%
[El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, N_el_domain, A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip] = Cavity_Mesh_Function(Receiver_Height,Receiver_Width,panel_n,panel_m,cover_n,cover_m,Plot_Mesh,Top_Lip,Bottom_Lip);


if Cal_View_Factors==1
    %%%Call F3D_30 View Factor (Panel to Panel view factors)
    [F_Panel_1_to_Panel] = F3D_30_Function(panel_m,panel_n,Receiver_Width,Receiver_Height,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel,Top_Lip,Bottom_Lip)

    %%%Call Monte Carlo Function - Panel 1 & 2 to floor 
    if cover_n==1 && cover_m==1
        [F_Panel_1_Floor  F_Panel_2_Floor F_Aper_Floor,F_Aper_Roof, F_Bottom_Lip_Floor, F_Top_Lip_Floor,  F_Bottom_Lip_Roof, F_Top_Lip_Roof] = Monte_Carlo_Panel_Floor_Pent_Function(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, Receiver_Width, Receiver_Height,Top_Lip, Bottom_Lip);
    else
        [F_Panel_1_Floor,F_Panel_2_Floor,F_Aper_Floor,F_Aper_Roof,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Bottom_Lip_Roof,F_Top_Lip_Roof] = MC_Panel_to_Floor_Function(Receiver_Height,Receiver_Width,N_rays,Top_Lip,Bottom_Lip,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel,El_Nod_Conn_Cover,Global_Nodal_Coor_Cover)
    end

    %%%Call Monte_Carlo Function Floor to Roof
    if cover_n==1 && cover_m==1
        [F_Roof_Floor] = MC_Roof_to_Floor_Pent_Function(Receiver_Height,Receiver_Width,N_rays,Global_Nodal_Coor_Cover);
    else
        [F_Roof_Floor] = MC_Roof_to_Floor_Function(Receiver_Width,Receiver_Height,cover_n,cover_m,N_rays,El_Nod_Conn_Cover,Global_Nodal_Coor_Cover)
    end

%     %%%Build View Factor Matrix
%     [VF_Matrix] = View_Factor_Matrix_Function_v3(N_el_domain,F_Panel_1_to_Panel,F_Panel_2_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename,panel_m);
% else
%      VF_Matrix=xlsread(Load_View_Factor_Matrix_File);
end

% VF_Matrix

% % %%%Vector of solar emissivities
% % for i=1:4*panel_m
% %    epsilon_sol(i)=e_act_sol; 
% % end
% % for i=1+panel_m*4:2+panel_m*4
% %     epsilon_sol(i)=e_pass_sol;
% % end
% % epsilon_sol(3+panel_m*4)=1;
% % 
% % %%%Vector of thermal emissivities
% % for i=1:4*panel_m
% %    epsilon_therm(i)=e_act_therm; 
% % end
% % for i=1+panel_m*4:2+panel_m*4
% %     epsilon_therm(i)=e_pass_therm;
% % end
% % epsilon_therm(3+panel_m*4)=1;
% % 
% % if Dist_Pattern==1
% %     Q_nodal_perc=ones(1,panel_m*4)*(1/panel_m);
% % elseif Dist_Pattern==2
% %     [i,j]=size(Q_nodal_perc);
% %     if j > panel_m || j < panel_m
% %         display("Error! Q_nodal_perc dimensions do not match number of nodes.")
% %     end
% %     temp=Q_nodal_perc;
% %     for i=2:4
% %         temp=[temp Q_nodal_perc]; 
% %     end
% %     Q_nodal_perc=temp;
% % else 
% %     [i,j]=size(Q_nodal_perc);
% %     if j > panel_m*4 || j < panel_m*4
% %         display("Error! Q_nodal_perc dimensions do not match number of nodes.")
% %     end    
% % end
% % 
% % 
% % %%%F Hat Matrix%%%
% % [F_hat_sol,F_hat_therm,epsilon_sol,epsilon_therm,rho_sol,rho_therm] = F_Hat_Function_v3(N_el_domain,VF_Matrix, epsilon_sol,epsilon_therm,Top_Lip,Bottom_Lip,e_pass_sol,e_pass_therm,panel_m);
% % 
% % %%%Energy Balance%%%
% % [T_max,T,Q_in_kW,Q_solar_loss_kW,Rad_Loss_kW,Conv_Loss_kW,HTF_Gain_kW,T_Nodal,h_bar_conv] = Energy_Balance_Function_v5(N_el_domain,epsilon_sol,epsilon_therm,rho_sol,rho_therm,F_hat_sol,F_hat_therm,T_HTF_in,T_HTF_out,T_infinity,A_Aperature,A_Panel,A_Roof,Receiver_Height,A_Top_Lip,A_Bottom_Lip,UA_HTF,h_Type,h_bar_conv,Q_in,panel_m,Top_Lip,Bottom_Lip,Q_nodal_perc);
% % 
% % A_Total=(A_Roof*2+A_Top_Lip*2);
% % Av_Temp_Passive=(T(5)*A_Roof/A_Total)+(T(6)*A_Roof/A_Total)+(T(7)*A_Top_Lip/A_Total)+(T(7)*A_Top_Lip/A_Total)
% % 
% % E_Balance_Cyl_Receiver_Function_v2(Receiver_Height,Receiver_Width,panel_m,Q_in,T_HTF_in,T_HTF_out,UA_HTF,e_act_sol,e_act_therm,T_infinity,h_bar_conv);

end

