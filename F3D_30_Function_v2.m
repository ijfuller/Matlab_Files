function [F_Panel] = F3D_30_Function_v2(panel_m,panel_n,Receiver_Width,Receiver_Height,El_Nod_Conn_Panel,Global_Nodal_Coord_Panel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%F3D-30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_el=panel_m*panel_n;
    interior_angle=135; %Angle between two adjacent panels
    panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel
    F_Panel=zeros(N_el,N_el,3);
    theta_set=[135,90,45, 67.5]; %Angles to panel 2, 3, 4, receiver aperature
    panel_int_offset=[0,cosd(theta_set(2)/2)*panel_width,(Receiver_Width/sind(45))*sind(67.5)-panel_width,0];
    for kk=1:4
        theta=theta_set(kk);
        for ii=1:N_el
            for jj=1:N_el
                Panel_1_Element=ii;
                Panel_2_Element=jj;

                if kk==4
                    z_1=0;
                    z_2=Receiver_Width;
                    eta_1=0;
                    eta_2=Receiver_Height;
                else
                    z_1=panel_width-Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_1_Element,2),1)+panel_int_offset(kk);
                    z_2=panel_width-Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_1_Element,1),1)+panel_int_offset(kk);
                    eta_1=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_1_Element,2),2);
                    eta_2=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_1_Element,3),2);
                end

                x_1=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_2_Element,1),1)+panel_int_offset(kk);
                x_2=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_2_Element,2),1)+panel_int_offset(kk);
                y_1=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_2_Element,2),2);
                y_2=Global_Nodal_Coord_Panel(El_Nod_Conn_Panel(Panel_2_Element,4),2);

                V(1,:)=[x_1 y_1 eta_1]; %Calulates G 8 times, V is the variable associated with the 8 runs
                V(2,:)=[x_1 y_1 eta_2];
                V(3,:)=[x_1 y_2 eta_1];
                V(4,:)=[x_1 y_2 eta_2];
                V(5,:)=[x_2 y_1 eta_1];
                V(6,:)=[x_2 y_1 eta_2];
                V(7,:)=[x_2 y_2 eta_1];
                V(8,:)=[x_2 y_2 eta_2];

                xi_1=z_1; %Formula in terms of xi but z is more understandable for coordinate systems
                xi_2=z_2;
                alpha=theta*pi/180; %Converts degree to radians for formula
                for i=1:8
                   x=V(i,1)+1e-8; %Adding a small value so that formula doesnt cause error
                   y=V(i,2)+1e-8;
                   eta=V(i,3); 
                   fun=@(xi) ((x-xi.*cos(alpha)).*cos(alpha)-xi.*sin(alpha).^2).*atan((eta-y)./sqrt(x.^2-2.*x.*xi.*cos(alpha)+xi.^2))./(sqrt(x.^2-2.*x.*xi.*cos(alpha)+xi.^2)*sin(alpha).^2)+cos(alpha).*(sqrt(xi.^2.*sin(alpha).^2+(eta-y).^2).*atan((x-xi.*cos(alpha))./sqrt(xi.^2*sin(alpha).^2+(eta-y).^2))-xi.*sin(alpha).*atan((x-xi.*cos(alpha))./sin(alpha)))./((eta-y).*sin(alpha).^2)+xi.*log((x.^2-2.*x.*xi.*cos(alpha)+xi.^2+(eta-y)^2)./(x.^2-2.*x.*xi.*cos(alpha)+xi.^2))./(2.*eta-2.*y);
                   G(i)=-(eta-y)*sin(alpha)^2/(2*pi)*integral(fun,xi_1,xi_2); %Calulates G for all 8 runs
                end
                %format long
                F_Panel(ii,jj,kk)=(-G(1)+G(2)+G(3)-G(4)+G(5)-G(6)-G(7)+G(8))/((x_2-x_1)*(y_2-y_1)); %Summation of 8 runs for view factor divided by area
            end
        end
    end
end

