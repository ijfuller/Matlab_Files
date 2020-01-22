function [F] = Monte_Carlo_Panel_Floor_Function(N_rays,El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover)

    figure(3)
    hold on
    grid on
    xlabel('X Axis') 
    ylabel('Y Axis') 
    zlabel('Z Axis') 
    xlim([0 20])
    ylim([0 12])
    zlim([0 14])

    for ii=1:length(El_Nod_Conn_Panel) %Panel loop
        %%%Panel Information
        rect_x=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,2),1)-Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,1),1); %Rectangular panel width
        rect_y=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,3),2)-Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,2),2); %Rectangular panel heigth
        offset_rect_y=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(ii,2),2); %Vertical offset

        plot3([0 rect_x],[0+offset_rect_y 0+offset_rect_y],[0 0],'k','Linewidth',2) %Bottom panel line
        plot3([rect_x rect_x],[0+offset_rect_y rect_y+offset_rect_y],[0 0],'k','Linewidth',2)%Right panel line
        plot3([rect_x 0],[rect_y+offset_rect_y rect_y+offset_rect_y],[0 0],'k','Linewidth',2)%Top panel line
        plot3([0 0],[rect_y+offset_rect_y 0+offset_rect_y],[0 0],'k','Linewidth',2)%Left panel line

        for jj=1:length(El_Nod_Conn_Cover) %Floor loop
            Element_MC=[ii jj]
            Pass_Fail_Set=zeros(1,4);
            %%%Floor Information
            floor_z=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),1) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),1)]+abs(Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(1,1),1));
            floor_y=zeros(1,4);
            floor_x=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),2) Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),2)];
            %offset_floor_z=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),1)+abs(Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(1,1),1));

            theta_rot=-22.5;
            floor_z_mod=floor_z.*cosd(theta_rot)+floor_x.*sind(theta_rot);
            floor_x_mod=-floor_z.*sind(theta_rot)+floor_x.*cosd(theta_rot);

            scatter3(floor_x_mod,floor_y,floor_z_mod,'filled','b')

            plot3([floor_x_mod(1) floor_x_mod(2)],[floor_y(1) floor_y(2)],[floor_z_mod(1) floor_z_mod(2)],'b','Linewidth',2) %Bottom floor line
            plot3([floor_x_mod(2) floor_x_mod(4)],[floor_y(2) floor_y(4)],[floor_z_mod(2) floor_z_mod(4)],'b','Linewidth',2) %Right floor line
            plot3([floor_x_mod(3) floor_x_mod(4)],[floor_y(3) floor_y(4)],[floor_z_mod(3) floor_z_mod(4)],'b','Linewidth',2) %Top floor line
            plot3([floor_x_mod(1) floor_x_mod(3)],[floor_y(1) floor_y(3)],[floor_z_mod(1) floor_z_mod(3)],'b','Linewidth',2) %Left floor line

            %%%Monte-Carlo Simulation
            ray_counter=0; %Initialize ray counter
            hit_counter=0; %Initialize hit counter
            x=[floor_z_mod(1), floor_z_mod(2), floor_z_mod(3), floor_z_mod(4)];
            y=[floor_x_mod(1), floor_x_mod(2), floor_x_mod(3), floor_x_mod(4)];
            m=[(y(2)-y(1))/(x(2)-x(1)), (y(4)-y(2))/(x(4)-x(2)), (y(3)-y(4))/(x(3)-x(4)), (y(1)-y(3))/(x(1)-x(3))]; %Bottom, Right Top Left
            y_int=[y(1)-m(1)*x(1), y(2)-m(2)*x(2), y(3)-m(3)*x(3), y(3)-m(4)*x(3)];
            while ray_counter<N_rays %Run loop until number of rays is reached
                ray_x=rand*rect_x; %Random x coordinate on panel
                ray_y=rand*rect_y+offset_rect_y; %Random y coordinate on panel

                Ptheta=rand; %Generate random number from 0 to 1
                theta=asin(sqrt(Ptheta)); %Calulate theta value
                Pphi=rand; %Generate random number from 0 to 1
                if Pphi<=0.5 || Pphi==1 %Conditions that result in automatic miss (ray would never hit floor element)
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

                if jj==1 | jj== 5 | jj== 9
                    [Pass_Fail_Set]=PIIQ_v1(m,y_int,x_random,y_random,x);
                elseif jj== 11 | jj== 12
                    [Pass_Fail_Set] = PIIQ_v2(m,y_int,x_random,y_random,y);
                else 
                    [Pass_Fail_Set]=PIIQ_v3(m,y_int,x_random,y_random);
                end
                    
%                 scatter3(ray_x,ray_y,0,'filled','k')
                if nnz(Pass_Fail_Set)==4 
                    hit_counter=hit_counter+1;
%                     plot3([ray_x,y_random],[ray_y,ry+ray_y],[0,x_random],'g')
%                     scatter3(y_random,ry+ray_y,x_random,'filled','g')
                else
%                     plot3([ray_x,y_random],[ray_y,ry+ray_y],[0,x_random],'r--')
%                     scatter3(y_random,ry+ray_y,x_random,'r')
                end
                
                ray_counter=ray_counter+1;
            end
        end
        F(ii,jj)=(hit_counter/ray_counter);
    end
end

display('------------------------------------------')
display('Monte Carlo Completed')
display(['Number of Rays per Simulation: ' num2str(N_rays)])
display('------------------------------------------')
end