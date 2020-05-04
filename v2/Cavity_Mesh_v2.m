clear all
clc

Receiver_Height=12;
Receiver_Width=14;
panel_n=1;
panel_m=1;
cover_n=6; %Elements in x direction. Less than 7
cover_m=3; %Elements in y direction. Less than 4
top_lip=1;
bottom_lip=1;

Plot_Mesh=1;

%%%%%%%%%%%%%Code Section%%%%%%%%%%%%%%%
interior_angle=135; %Angle between two adjacent panels
panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel
El_Panel=panel_n*panel_m; %Number of elements for a single panel
El_Cover=cover_n*cover_m; %Number of elements for the roof/floor
N_el_domain=4*(panel_n*panel_m)+2*(cover_n*cover_m)+1; %Number of elements in domain

% A_Roof=(1+sqrt(2))*panel_width^2; %Area of roof elements
% A_Panel=(panel_width*Receiver_Height); %Area of panels
% A_Aperature=Receiver_Height*Receiver_Width; %Area of aperature


%%%%%For panels%%%%%%%%
% panel_x=[0:(panel_width/panel_n):panel_width]; %Panel nodal locations local in x axis
% panel_y=[0:(Receiver_Height/panel_m):Receiver_Height]; %Panel nodal locations local in y axis
% 
% %%%Global Coordinates Panels%%%
% Panels=4; %Number of panels
% Global_Nodal_Coord_Panel=zeros((panel_n+1)*(panel_m+1),3); %Initialize global nodal coordinates x, y, z for one panel.
% x_coord=[0:(panel_width/panel_n):panel_width]; %X nodal location of panel
% y_coord=[0:(Receiver_Height/panel_m):Receiver_Height]; %Y nodal location of panels
% 
% k=1; %Start filling in row one of global panel coordinates
% for i=1:panel_m+1 %Loop over y coordinates
%    for j=1:panel_n+1 %Loop over x coordinates
%       Global_Nodal_Coord_Panel(k,:,1)=[x_coord(j),y_coord(i),0]; %Add coordinates to each row, x,y,z
%       k=k+1; %Move to next row
%    end
% end
% 
% N_nod_x=panel_n+1; %Number of nodes in x direction of panel
% N_nod_y=panel_m+1; %Number of nodes in y direction of panel
% Nod_per_el=4; %Nodes per quadrilateral element
% N_el=panel_n*panel_m; %Number of elements in domain for quadralateral
% El_Nod_Conn_Panel=zeros(N_el,Nod_per_el); %Initializing the node connectivity matrix
% 
% i=1;%Start at row one of connectivity matrix
% for j=1:panel_m
%     for k=1:panel_n
%         El_Nod_Conn_Panel(i,:)=[k+(j-1)*N_nod_x, k+1+(j-1)*N_nod_x,k+N_nod_x+(j-1)*N_nod_x, k+N_nod_x+1+(j-1)*N_nod_x]; %Node connectivity for every element in one panel
%         i=i+1; %Go to next row
%     end
% end
%%%%%For panels end%%%%%%%%

%%%%Roof and Floor of Cavity Start%%%%
panel_angle=[-(180-(interior_angle/2)-90)*3 -(180-(interior_angle/2)-90) (180-(interior_angle/2)-90) (180-(interior_angle/2)-90)*3]; %Calulates the angle rotation from centerline for panels
panel_radial_dist=cosd(panel_angle(2))*(Receiver_Width/2); %Distance between the cavity origin (middle of aperature) and middle of panel
roof_x=[(-Receiver_Width/2) (-Receiver_Width/2)*cosd(45) 0 (Receiver_Width/2)*cosd(45) (Receiver_Width/2)]; %X coordinate of roof
roof_y=[0 (Receiver_Width/2)*cosd(45) (Receiver_Width/2) (Receiver_Width/2)*cosd(45) 0]; %Y coordinate of roof

 
if cover_m==1 && cover_n==1 %If roof is a single element
    Global_Nodal_Coor_Cover=[roof_x;roof_y]'; %X and Y coordinates of each 5 points
    El_Nod_Conn_Cover=[1 5 4 3 2]; %Set element connectivity of cover
else
    %%%Top and Bottom Mesh%%%
    x_div_loc=[-Receiver_Width/2:(Receiver_Width/(cover_n)):Receiver_Width/2]; %Location for x divisions
    y_div_loc=[0:((Receiver_Width/2)/cover_m):(Receiver_Width/2)]; %Location of y divisions

    %which panel does each x division line intersect
    for i=1:(length(x_div_loc)-2)
        if roof_x(1)<=x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(2)
            Loc(i)=1; %Panel 1
        end
        if roof_x(2)<x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(3)
            Loc(i)=2; %Panel 2
        end
        if roof_x(3)<x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(4)
            Loc(i)=3; %Panel 3
        end
        if roof_x(4)<x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(5)
            Loc(i)=4; %Panel 4
        end
    end

    %%%Equations for panels lines. Slopes and y intercepts.
    m(1)=(roof_y(1)-roof_y(2))/(roof_x(1)-roof_x(2));
    y_int(1)=roof_y(1)-m(1)*roof_x(1);
    m(2)=(roof_y(2)-roof_y(3))/(roof_x(2)-roof_x(3));
    y_int(2)=roof_y(2)-m(2)*roof_x(2);
    m(3)=(roof_y(3)-roof_y(4))/(roof_x(3)-roof_x(4));
    y_int(3)=roof_y(3)-m(3)*roof_x(3);
    m(4)=(roof_y(4)-roof_y(5))/(roof_x(4)-roof_x(5));
    y_int(4)=roof_y(4)-m(4)*roof_x(4);

    
    
    if Plot_Mesh==1 %%%Plot roof external points
        figure(2)
        hold on
        scatter(roof_x,roof_y,'k','filled')
        plot(roof_x,roof_y,'k','LineWidth',2)
        xlim([-Receiver_Width/2 Receiver_Width/2])
        ylim([0 Receiver_Width])
        %%%Get four roof external points and plot them end
        title('Roof & Floor Discretization','FontSize', 28)
        xlabel('X position [m]','FontSize', 24) 
        ylabel('Y Position [m]','FontSize', 24) 
    end
    
    
    %%%New Code Addition%%%
    if cover_n>1
       for i=1:(cover_n-1)
           X_Loc=x_div_loc(i+1);
           Y_Loc=m(Loc(i))*x_div_loc(i+1)+y_int(Loc(i));
           plot([X_Loc X_Loc],[0 Y_Loc],'k','LineWidth',1)
       end
    end

%     %which panel does each y division line intersect
%     Loc=0; %Reset panel intersection value to zero
%     if 1<cover_m
%         for i=1:(length(y_div_loc)-2)
%             if roof_y(1)<=y_div_loc(i+1) && y_div_loc(i+1)<=roof_y(2)
%                 Loc(i)=1; %Panel 1
%             end
%             if roof_y(2)<y_div_loc(i+1) && y_div_loc(i+1)<=roof_y(3)
%                 Loc(i)=2; %Panel 2
%             end
%         end
%         for i=1:(cover_m-1)
%            Y_loc=y_div_loc(i+1);
%            X_loc=(Y_loc-y_int(Loc(i)))/m(Loc(i));
%            plot([X_loc -X_loc],[Y_loc Y_loc],'k','LineWidth',1)
%         end
%     end

if cover_n>1
        for i=1:(length(x_div_loc)-2)
            height(i)=m(Loc(i))*x_div_loc(i+1)+y_int(Loc(i));
        end

        for i=1:length(height)
            y_div_loc_line(i,:)=[0:(height(i)/cover_m):height(i)];
        end
    end
    y_cover_coord=zeros(cover_n*2,cover_m+1);
    j=1;
    for i=1:cover_n/2
        y_cover_coord(j,:)=y_div_loc_line(i,:);
        y_cover_coord(j+1,:)=y_div_loc_line(i,:);
        j=j+2;
    end
    if mod(cover_n,2)==1
        y_cover_coord(j,:)=y_div_loc_line(floor(cover_n/2),:);
    end
    y_cover_coord(1,cover_n)=(Receiver_Width/2)*cosd(45);
    y_cover_coord=y_cover_coord+flipud(y_cover_coord);

    x_cover_coord=zeros(cover_n*2,cover_m+1);
    j=1;
    for i=1:cover_n/2
        x_cover_coord(j,:)=ones(1,cover_m+1)*x_div_loc(i);
        x_cover_coord(j+1,:)=ones(1,cover_m+1)*x_div_loc(i+1);
        j=j+2;
    end
    if mod(cover_n,2)==1
        x_cover_coord(j,:)=x_div_loc(floor(cover_n/2),:);
    end
    x_cover_coord=x_cover_coord+flipud(x_cover_coord)*-1;
    for i=1:cover_m+1
       x_cover_coord(1,i)=(y_cover_coord(1,i)-y_int(1))/m(1);
       x_cover_coord(cover_n*2,i)=(y_cover_coord(cover_n*2,i)-y_int(4))/m(4);
    end

        El_Nod_Conn_Cover=zeros(El_Cover,4);
    k=1;
    for i=1:2:El_Cover*2
        El_Nod_Conn_Cover(k,1)=i;
        El_Nod_Conn_Cover(k,2)=i+1;
        El_Nod_Conn_Cover(k,3)=i+(2*cover_n);
        El_Nod_Conn_Cover(k,4)=i+1+(2*cover_n);
        k=k+1;
    end
    %%%End of roof descritization

    %%%%Global Coordinates%%%%
    [i j]=size(x_cover_coord);
    Global_Nodal_Coor_Cover=reshape(x_cover_coord,[i*j,1]);
    Global_Nodal_Coor_Cover(:,2)=reshape(y_cover_coord,[i*j,1]);
    
    %%%New Code Addition End%%%
end









%%%Plotting%%%
% if Plot_Mesh==1 %%%Plot panel elements
%     figure(1)
%     hold on
%     %%%plots nodal points of panel
%     for i=1:length(panel_x)
%         for j=1:length(panel_y)
%             scatter(panel_x(i),panel_y(j),'k','filled')
%         end
%     end
%     xlim([0 panel_width])
%     ylim([0 Receiver_Height])
% 
%     %%%plot lines of panel
%     for i=1:length(panel_x)
%         line_x=ones(1,length(panel_y))*panel_x(i);
%         plot(line_x,panel_y,'k')
%     end
%     for i=1:length(panel_y)
%         line_y=ones(1,length(panel_x))*panel_y(i);
%         plot(panel_x,line_y,'k')
%     end
%     title('Panel Discretization','FontSize', 28)
%     xlabel('X position [m]','FontSize', 24) 
%     ylabel('Y Position [m]','FontSize', 24) 
% end


 
% if Plot_Mesh==1
%     %%%Plot Roof Descritization
%     for i=1:length(Loc)
%         plot([x_div_loc(i+1),x_div_loc(i+1)],[0,height(i)],'k')
%     end
%     for i=1:2:cover_n*2
%         for j=2:cover_m
%             plot([x_cover_coord(i,j),x_cover_coord(i+1,j)],[y_cover_coord(i,j),y_cover_coord(i+1,j)],'k')
%         end
%     end
% end


% display('------------------------------------------')
% display('Cavity Mesh Generation Completed')
% display(['Panel Elements:' num2str(El_Panel) ' '])
% display(['Roof/Floor Elements:' num2str(El_Cover) ' '])
% display(['Total Number of Elements:' num2str(N_el_domain) ' '])
% display('------------------------------------------')