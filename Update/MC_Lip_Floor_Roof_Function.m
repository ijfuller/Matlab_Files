function [F_Bottom_Lip_Floor,F_Top_Lip_Floor,F_Bottom_Lip_Roof,F_Top_Lip_Roof] = MC_Lip_Floor_Roof_Function(Receiver_Height,Receiver_Width,Top_Lip,Bottom_Lip,N_rays,El_Nod_Conn_Panel,El_Nod_Conn_Cover,Global_Nodal_Coord_Panel,Global_Nodal_Coor_Cover,cover_n,cover_m)
 N_Lip_Runs=2;

    if Bottom_Lip > 0
        N_Lip_Runs=N_Lip_Runs+2;
    end

    plot_MC_Ray=0; %Disable=0, Enable=1. This will plot the rays/hits/misses. For code modification only.

    plot_panel=1; %Which panel to plot. Outer panel (1), inner panel (2), aperature-floor (3), top lip-floor (5), top lip-roof (6), bottom lip-floor (7), bottom lip-roof (8)

    F_Bottom_Lip_Floor=0; %Initialize
    F_Top_Lip_Floor=0; %Initialize
    F_Bottom_Lip_Roof=0; %Initialize
    F_Top_Lip_Roof=0; %Initialize

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

    for ii=1:N_Lip_Runs
        offset_rect_x=0;

        if ii==1
                %Top Lip to floor
                    rect_x=Receiver_Width;
                    rect_y=Top_Lip;
                    offset_rect_y=Bottom_Lip+(Receiver_Height-Bottom_Lip-Top_Lip);
        end
        if ii==2
                %Top Lip to roof
                    rect_x=Receiver_Width;
                    rect_y=Top_Lip;
                    offset_rect_y=0;
        end 
        if ii==3; %Bottom Lip to floor
            rect_x=Receiver_Width;
            rect_y=Bottom_Lip;
            offset_rect_y=0;
        end
        if ii==4; %Bottom Lip to roof
            rect_x=Receiver_Width;
            rect_y=Bottom_Lip;
            offset_rect_y=+Top_Lip+(Receiver_Height-Bottom_Lip-Top_Lip);
        end

            
            if plot_MC_Ray==1 && panel==plot_panel;
                plot3([0+offset_rect_x rect_x+offset_rect_x],[0+offset_rect_y 0+offset_rect_y],[0 0],'k','Linewidth',2) %Bottom panel line
                plot3([rect_x+offset_rect_x rect_x+offset_rect_x],[0+offset_rect_y rect_y+offset_rect_y],[0 0],'k','Linewidth',2)%Right panel line
                plot3([rect_x+offset_rect_x 0+offset_rect_x],[rect_y+offset_rect_y rect_y+offset_rect_y],[0 0],'k','Linewidth',2)%Top panel line
                plot3([0+offset_rect_x 0+offset_rect_x],[rect_y+offset_rect_y 0+offset_rect_y],[0 0],'k','Linewidth',2)%Left panel line
            end

            for jj=1:(cover_n*cover_m)
                Pass_Fail_Set=zeros(1,4);
                %%%Floor Information
                floor_z=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),1),Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),1),Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),1),Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),1)]+abs(Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(1,1),1));
                floor_y=zeros(1,4);
                floor_x=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),2),Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),2),Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),2),Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),2)]+abs(Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(1,1),2));

                theta_rot=90; %Rotation for opening to floor
                floor_z_mod=floor_z.*cosd(theta_rot)+floor_x.*sind(theta_rot);
                floor_x_mod=-floor_z.*sind(theta_rot)+floor_x.*cosd(theta_rot)+Receiver_Width;
                
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
                x=[floor_z_mod(1), floor_z_mod(2), floor_z_mod(3), floor_z_mod(4)];
                y=[floor_x_mod(1), floor_x_mod(2), floor_x_mod(3), floor_x_mod(4)];
                m=[(y(2)-y(1))/(x(2)-x(1)), (y(3)-y(2))/(x(3)-x(2)), (y(4)-y(3))/(x(4)-x(3)), (y(1)-y(4))/(x(1)-x(4))]; %Bottom, Right Top Left
                y_int=[y(1)-m(1)*x(1), y(2)-m(2)*x(2), y(3)-m(3)*x(3), y(4)-m(4)*x(4)];
                while ray_counter<N_rays %Run loop until number of rays is reached
                    ray_x=rand*rect_x+offset_rect_x; %Random x coordinate on panel
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


                    [Pass_Fail_Set] = PIIQ_90_Function(m,y_int,x_random,y_random,x,y,El_Nod_Conn_Cover,Global_Nodal_Coor_Cover,jj);
         

                    if plot_MC_Ray==1 && panel==plot_panel; %Plot
                        scatter3(ray_x,ray_y,0,'filled','k')
                    end
                    
                    if nnz(Pass_Fail_Set)==4; 
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
    if ii==1
        F_Top_Lip_Floor=F(1,:); 
    end
    if ii==2
        F_Top_Lip_Roof=F(2,:); 
    end
    if ii==3
        F_Bottom_Lip_Floor=F(3,:); 
    end
    if ii==4
        F_Bottom_Lip_Roof=F(4,:); 
    end   
    display('-----------------------------------------------')
    display('Monte Carlo Lip(s)-Floor/Roof Completed')
    display(['Number of Rays per Simulation: ' num2str(N_rays)])
    display('-----------------------------------------------')
end