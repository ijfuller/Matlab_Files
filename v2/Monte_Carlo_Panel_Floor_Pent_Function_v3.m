function [F_Panel_1_Floor  F_Panel_2_Floor F_Aper_Floor, F_Bottom_Lip_Floor, F_Top_Lip_Floor] = Monte_Carlo_Panel_Floor_Pent_Function_v3(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, Receiver_Width, Receiver_Height, Top_Lip, Bottom_Lip)

plot_MC_Ray=1; %Disable=0, Enable=1. This will plot the rays/hits/misses. For code modification only.

plot_panel=5; %Which panel to plot. Outer panel (1), inner panel (2), aperature (3), top lip-floor (4), bottom lip-floor (5)

if plot_MC_Ray==1;
    figure(3)
    hold on
    grid on
    xlabel('X Axis') 
    ylabel('Y Axis') 
    zlabel('Z Axis') 
    xlim([-5 20])
    ylim([0 12])
    zlim([0 14])
end
    
    for panel=1:5 %Think this is panel 1, panel 2, aperature
        [h w]=size(El_Nod_Conn_Panel);
        for ii=1:h %Loop over number of elements/panel
            %%%Panel Information
            rect_x=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,2),1)-Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,1),1); %Rectangular panel width
            rect_y=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,3),2)-Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,2),2); %Rectangular panel heigth
            offset_rect_y=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,2),2); %Vertical offset
            if panel==3; %Aperature
                rect_x=Receiver_Width;
                rect_y=Receiver_Height-Top_Lip-Bottom_Lip;
                offset_rect_y=Bottom_Lip;
            end
            if panel==4; %Top Lip
                rect_x=Receiver_Width;
                rect_y=Top_Lip;
                offset_rect_y=Bottom_Lip+(Receiver_Height-Bottom_Lip-Top_Lip);
            end
            if panel==5; %Bottom Lip
                rect_x=Receiver_Width;
                rect_y=Bottom_Lip;
                offset_rect_y=0;
            end
            
            if plot_MC_Ray==1 && panel==plot_panel;
                plot3([0 rect_x],[0+offset_rect_y 0+offset_rect_y],[0 0],'k','Linewidth',2) %Bottom panel line
                plot3([rect_x rect_x],[0+offset_rect_y rect_y+offset_rect_y],[0 0],'k','Linewidth',2)%Right panel line
                plot3([rect_x 0],[rect_y+offset_rect_y rect_y+offset_rect_y],[0 0],'k','Linewidth',2)%Top panel line
                plot3([0 0],[rect_y+offset_rect_y 0+offset_rect_y],[0 0],'k','Linewidth',2)%Left panel line
            end

            for jj=1
                Pass_Fail_Set=zeros(1,5);
                %%%Floor Information
                floor_z=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,5),1)]+abs(Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(1,1),1));
                floor_y=zeros(1,5);
                floor_x=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,5),2)];

            if panel==1;
                theta_rot=-22.5; %For panel 1
                floor_z_mod=floor_z.*cosd(theta_rot)+floor_x.*sind(theta_rot);
                floor_x_mod=-floor_z.*sind(theta_rot)+floor_x.*cosd(theta_rot);
            elseif panel==3 || panel==4 || panel==5;
                theta_rot=90; %Rotation for opening to floor
                floor_z_mod=floor_z.*cosd(theta_rot)+floor_x.*sind(theta_rot);
                floor_x_mod=-floor_z.*sind(theta_rot)+floor_x.*cosd(theta_rot)+Receiver_Width;
            else
                theta_rot=-67.5; %%%Panel 2
                floor_z_mod=floor_z.*cosd(theta_rot)+floor_x.*sind(theta_rot);
                floor_x_mod=-floor_z.*sind(theta_rot)+floor_x.*cosd(theta_rot);
                floor_z_mod=floor_z_mod+abs(floor_z_mod(5));
                floor_x_mod=floor_x_mod-abs(floor_x_mod(5));
            end

            if plot_MC_Ray==1 && panel==plot_panel;
                scatter3(floor_x_mod,floor_y,floor_z_mod,'filled','b')
                plot3([floor_x_mod(1) floor_x_mod(2)],[floor_y(1) floor_y(2)],[floor_z_mod(1) floor_z_mod(2)],'b','Linewidth',2) 
                plot3([floor_x_mod(2) floor_x_mod(3)],[floor_y(2) floor_y(3)],[floor_z_mod(2) floor_z_mod(3)],'b','Linewidth',2) 
                plot3([floor_x_mod(3) floor_x_mod(4)],[floor_y(3) floor_y(4)],[floor_z_mod(3) floor_z_mod(4)],'b','Linewidth',2) 
                plot3([floor_x_mod(4) floor_x_mod(5)],[floor_y(4) floor_y(5)],[floor_z_mod(4) floor_z_mod(5)],'b','Linewidth',2) 
                plot3([floor_x_mod(1) floor_x_mod(5)],[floor_y(1) floor_y(5)],[floor_z_mod(1) floor_z_mod(5)],'b','Linewidth',2)
            end

                %%%Monte-Carlo Simulation
                ray_counter=0; %Initialize ray counter
                hit_counter=0; %Initialize hit counter
                x=[floor_z_mod(1), floor_z_mod(2), floor_z_mod(3), floor_z_mod(4) floor_z_mod(5)];
                y=[floor_x_mod(1), floor_x_mod(2), floor_x_mod(3), floor_x_mod(4) floor_x_mod(5)];
                m=[(y(2)-y(1))/(x(2)-x(1)), (y(3)-y(2))/(x(3)-x(2)), (y(4)-y(3))/(x(4)-x(3)), (y(5)-y(4))/(x(5)-x(4)), (y(1)-y(5))/(x(1)-x(5))]; %Bottom, Right Top Left
                y_int=[y(1)-m(1)*x(1), y(2)-m(2)*x(2), y(3)-m(3)*x(3), y(4)-m(4)*x(4) , y(5)-m(5)*x(5)];
                while ray_counter<N_rays %Run loop until number of rays is reached
                    ray_x=rand*rect_x; %Random x coordinate on panel
                    ray_y=rand*rect_y+offset_rect_y; %Random y coordinate on panel

                    Ptheta=rand; %Generate random number from 0 to 1
                    theta=asin(sqrt(Ptheta)); %Calulate theta value
                    Pphi=rand; %Generate random number from 0 to 1
                    if Pphi<=0.5 || Pphi==1; %Conditions that result in automatic miss (ray would never hit floor element)
                        ray_counter=ray_counter+1;
                    else 
                    phi=Pphi*2*pi; %Calculate phi value

                    rx=cos(phi)*sin(theta); %Unit x component
                    ry=sin(phi)*sin(theta); %Unit y component
                    rz=cos(theta); %Unit z component
                    scal=abs(ry)/ray_y; %Scalar value for unit vector 
                    rx=rx/scal; %Intersection x component
                    ry=ry/scal; %Intersection y component
                    rz=rz/scal; %Intersection z component

                    mul=abs(ry)/ray_y;
                    rx=rx/mul;
                    ry=ry/mul;
                    rz=rz/mul;

                    x_random=rz; %point on which the ray strikes
                    y_random=rx+ray_x;
                    
                    if panel==1;
                        [Pass_Fail_Set] = PI_Pent_Function(m,y_int,x_random,y_random, y);
                    elseif panel==3 || panel==4 || panel==5
                       [Pass_Fail_Set] = PI_Pent_Floor_Function(m,y_int,x_random,y_random, x, y);
                    else
                        [Pass_Fail_Set] = PI_Pent_Pan_2_Function(m,y_int,x_random,y_random, x, y);
                    end

                    if plot_MC_Ray==1 && panel==plot_panel; %Plot
                        scatter3(ray_x,ray_y,0,'filled','k')
                    end
                    
                    if nnz(Pass_Fail_Set)==5; 
                        hit_counter=hit_counter+1;
                        if plot_MC_Ray==1 && panel==plot_panel;%Plot
                            plot3([ray_x,y_random],[ray_y,ry+ray_y],[0,x_random],'g')
                            scatter3(y_random,ry+ray_y,x_random,'filled','g')
                        end
                    else
                        if plot_MC_Ray==1 && panel==plot_panel;%Plot
                            plot3([ray_x,y_random],[ray_y,ry+ray_y],[0,x_random],'r--')
                            scatter3(y_random,ry+ray_y,x_random,'r')
                        end
                    end

                    ray_counter=ray_counter+1;
                end
            end
            F(ii,jj)=(hit_counter/ray_counter);
        end
        end
        if panel==1;
            F_Panel_1_Floor=F;
        elseif panel==3
            F_Aper_Floor=F;
        elseif panel==4
            F_Top_Lip_Floor=F;
        elseif panel==5
            F_Bottom_Lip_Floor=F;
        else panel==2;
            F_Panel_2_Floor=F;
        end

    end
    display('------------------------------------------')
    display('Monte Carlo Completed')
    display(['Number of Rays per Simulation: ' num2str(N_rays)])
    display('------------------------------------------')
end