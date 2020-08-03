function [El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, N_el_domain, A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip] = Cavity_Mesh_Function(Receiver_Height,Receiver_Width,panel_n,panel_m,cover_n,cover_m,Plot_Mesh,Top_Lip,Bottom_Lip)

%%%%% Function Meshes the panels and roof/floor %%%%%
interior_angle=135; %Angle between two adjacent panels in degrees.
panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel in meters
El_Panel=panel_n*panel_m; %Number of elements for a single panel
El_Cover=cover_n*cover_m; %Number of elements for the roof/floor
N_el_domain=4*(panel_n*panel_m)+2*(cover_n*cover_m)+1; %Number of elements in domain

if Top_Lip > 0
    N_el_domain=N_el_domain+1; %Add one element if the top lip is active.
end
if Bottom_Lip > 0
    N_el_domain=N_el_domain+1; %Adds one element if the bottom lip is active.
end

A_Roof=(1+sqrt(2))*panel_width^2; %Area of the roof/floor
A_Panel=panel_width*Receiver_Height; %Area of a single panel
A_Top_Lip=Top_Lip*Receiver_Width; %Area of Top Lip
A_Bottom_Lip=Bottom_Lip*Receiver_Width; %Area of bottom lip
A_Aperature=(Receiver_Height*Receiver_Width)-A_Top_Lip-A_Bottom_Lip; %Area of aperature

%%%%% Panel Mesh Section %%%%%
Panels=4; %Number of panels
N_nod_x=panel_n+1; %Number of nodes in x direction
N_nod_y=panel_m+1; %Number of nodes in y direction
Nod_per_el=4; %Nodes per quadrilateral element
N_el=panel_n*panel_m; %Number of elements in domain for quadralateral
panel_x=[0:(panel_width/panel_n):panel_width]; %Panel nodal locations local in x axis
panel_y=[0:(Receiver_Height/panel_m):Receiver_Height]; %Panel nodal locations local in y axis
Global_Nodal_Coord_Panel=zeros(panel_n*panel_m,3); %Intializes nodal coordinates (X,Y,Z)
x_coord=[0:(panel_width/panel_n):panel_width]; %X coordinate location of nodes
y_coord=[0:(Receiver_Height/panel_m):Receiver_Height]; %Y coordinate location of nodes
k=1; %Start at row 1
for i=1:N_nod_y %Loop over vertical nodes
   for j=1:N_nod_x %Loop over horizontal nodes
      Global_Nodal_Coord_Panel(k,:,1)=[x_coord(j),y_coord(i),0]; %Fills in coordinates for every node of a single panel
      k=k+1; %Go to next row
   end
end
El_Nod_Conn_Panel=zeros(N_el,Nod_per_el); %Initializing the node connectivity matrix
i=1; %Start at row 1
for j=1:panel_m %Loop over vertical elements
    for k=1:panel_n %Loop over horizontal elements
        El_Nod_Conn_Panel(i,:)=[k+(j-1)*N_nod_x, k+1+(j-1)*N_nod_x,k+N_nod_x+(j-1)*N_nod_x, k+N_nod_x+1+(j-1)*N_nod_x]; %Node connectivity for every element
        i=i+1; %Go to next row
    end
end
%%%Plot Panel Mesh%%%
if Plot_Mesh==1
    figure(1) 
    hold on
    %%%plots nodal points of panel
    for i=1:length(panel_x)
        for j=1:length(panel_y)
            scatter(panel_x(i),panel_y(j),'k','filled') %Fill in panel nodes
        end
    end
    xlim([0 panel_width]) %Set x limit of figure
    ylim([0 Receiver_Height]) %Set y limit of figure

    %%%plot lines of panel
    for i=1:length(panel_x)
        line_x=ones(1,length(panel_y))*panel_x(i); %Horizontal lines
        plot(line_x,panel_y,'k') %Plot horizontal lines
    end
    for i=1:length(panel_y)
        line_y=ones(1,length(panel_x))*panel_y(i); %Vertical lines
        plot(panel_x,line_y,'k') %Plot vertical lines
    end
    title('Panel Discretization','FontSize', 28)
    xlabel('X position [m]','FontSize', 24) 
    ylabel('Y Position [m]','FontSize', 24) 
    hold off
end
%%%%% End of Panel Mesh Section %%%%%

%%%Cavity Top & Bottom Cover
%%%Get four roof external points
panel_angle=[-(180-(interior_angle/2)-90)*3 -(180-(interior_angle/2)-90) (180-(interior_angle/2)-90) (180-(interior_angle/2)-90)*3]; %Calulates the angle rotation from centerline for panels
panel_radial_dist=cosd(panel_angle(2))*(Receiver_Width/2); %Distance between the cavity origin and middle of panel
roof_x=[(-Receiver_Width/2) (-Receiver_Width/2)*cosd(45) 0 (Receiver_Width/2)*cosd(45) (Receiver_Width/2)]; %X coordinates of external power nodes
roof_y=[0 (Receiver_Width/2)*cosd(45) (Receiver_Width/2) (Receiver_Width/2)*cosd(45) 0]; %Y coordinates of external power nodes

m(1)=(roof_y(1)-roof_y(2))/(roof_x(1)-roof_x(2)); %Slope of line 1
y_int(1)=roof_y(1)-m(1)*roof_x(1); %Y intercept of line 1
m(2)=(roof_y(2)-roof_y(3))/(roof_x(2)-roof_x(3)); %Slope of line 2
y_int(2)=roof_y(2)-m(2)*roof_x(2);%Y intercept of line 2
m(3)=(roof_y(3)-roof_y(4))/(roof_x(3)-roof_x(4)); %Slope of line 3
y_int(3)=roof_y(3)-m(3)*roof_x(3);%Y intercept of line 3
m(4)=(roof_y(4)-roof_y(5))/(roof_x(4)-roof_x(5)); %Slope of line 4
y_int(4)=roof_y(4)-m(4)*roof_x(4);%Y intercept of line 4

%%%Element and nodal coordinates for various meshes%%%
if cover_n==1 && cover_m==1
   Global_Nodal_Coor_Cover=[roof_x;roof_y]';
   El_Nod_Conn_Cover=[1 5 4 3 2];
end

if cover_n==2 && cover_m==1
    Global_Nodal_Coor_Cover=[(-Receiver_Width/2),0; 0,0; (Receiver_Width/2),0; (-Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45); 0,(Receiver_Width/2); (Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45)];
    El_Nod_Conn_Cover=[1 2 5 4; 2 3 6 5];
end

if cover_n==4 && cover_m==1
    Global_Nodal_Coor_Cover=[(-Receiver_Width/2),0; (-Receiver_Width/4),0; 0,0; (Receiver_Width/4),0; (Receiver_Width/2),0; (-Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45); -(Receiver_Width/4),(m(2)*-(Receiver_Width/4)+y_int(2)); 0,(Receiver_Width/2); (Receiver_Width/4),(m(3)*(Receiver_Width/4)+y_int(3)); (Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45)];
    El_Nod_Conn_Cover=[1 2 7 6; 2 3 8 7; 3 4 9 8; 4 5 10 9];
end

if cover_n==2 && cover_m==2
    Global_Nodal_Coor_Cover=[(-Receiver_Width/2),0; 0,0; (Receiver_Width/2),0; -(m(2)*-(Receiver_Width/4)+y_int(2)),(Receiver_Width/4); 0,(Receiver_Width/4); (m(2)*-(Receiver_Width/4)+y_int(2)),(Receiver_Width/4); (-Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45); 0,(Receiver_Width/2); (Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45)];
    El_Nod_Conn_Cover=[1 2 5 4; 2 3 6 5; 4 5 8 7; 5 6 9 8];
end

if cover_n==2 && cover_m==3
    Global_Nodal_Coor_Cover=[(-Receiver_Width/2),0; 0,0; (Receiver_Width/2),0; ((Receiver_Width/6)-y_int(1))/m(1),(Receiver_Width/6); 0,(Receiver_Width/6); ((Receiver_Width/6)-y_int(4))/m(4),(Receiver_Width/6); ((Receiver_Width/3)-y_int(1))/m(1),(Receiver_Width/3); 0,(Receiver_Width/3); ((Receiver_Width/3)-y_int(4))/m(4),(Receiver_Width/3); -(Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45); 0,(Receiver_Width/2); (Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45)];
    El_Nod_Conn_Cover=[1 2 5 4; 2 3 6 5; 4 5 8 7; 5 6 9 8; 7 8 11 10; 8 9 12 11];
end

if cover_n==4 && cover_m==3
    Global_Nodal_Coor_Cover=[(-Receiver_Width/2),0; (-Receiver_Width/4),0; 0,0; (Receiver_Width/4),0; (Receiver_Width/2),0;(((m(2)*(-Receiver_Width/4)+y_int(2))/3)-y_int(1))/m(1),(m(2)*((-Receiver_Width/4))+y_int(2))/3; (-Receiver_Width/4),(m(2)*((-Receiver_Width/4))+y_int(2))/3; (-Receiver_Width/4),(Receiver_Width/2)/3;0,(Receiver_Width/2)/3; (Receiver_Width/4),(Receiver_Width/2)/3;(Receiver_Width/4),(m(3)*((Receiver_Width/4))+y_int(3))/3; (((m(3)*(Receiver_Width/4)+y_int(3))/3)-y_int(4))/m(4),(m(3)*((Receiver_Width/4))+y_int(3))/3;(((m(2)*(-Receiver_Width/4)+y_int(2))*2/3)-y_int(1))/m(1),(m(2)*((-Receiver_Width/4))+y_int(2))*2/3;(-Receiver_Width/4),(m(2)*((-Receiver_Width/4))+y_int(2))*2/3;(-Receiver_Width/4),((Receiver_Width/2)*2)/3;0,((Receiver_Width/2)*2)/3;(Receiver_Width/4),((Receiver_Width/2)*2)/3;(Receiver_Width/4),(m(2)*((-Receiver_Width/4))+y_int(2))*2/3;(((m(3)*(Receiver_Width/4)+y_int(3))*2/3)-y_int(4))/m(4),(m(3)*((Receiver_Width/4))+y_int(3))*2/3;(-Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45);(-Receiver_Width/4),m(2)*(-Receiver_Width/4)+y_int(2);0, (Receiver_Width/2);(Receiver_Width/4),m(3)*(Receiver_Width/4)+y_int(3);(Receiver_Width/2)*cosd(45),(Receiver_Width/2)*cosd(45)];
    El_Nod_Conn_Cover=[1 2 7 6; 2 3 9 8; 3 4 10 9; 4 5 12 11; 6 7 14 13; 8 9 16 15; 9 10 17 16; 11 12 19 18; 13 14 21 20; 15 16 22 21; 16 17 23 22; 18 19 24 23];
end

%%%Plot Mesh for floor/roof%%%
if Plot_Mesh==1
    figure(2)
    hold on
    [Len,Wid]=size(El_Nod_Conn_Cover);
    if cover_n==1 && cover_m==1
        plot([roof_x(1),roof_x(2)],[roof_y(1),roof_y(2)],'k','Linewidth',2)    
        plot([roof_x(2),roof_x(3)],[roof_y(2),roof_y(3)],'k','Linewidth',2)  
        plot([roof_x(3),roof_x(4)],[roof_y(3),roof_y(4)],'k','Linewidth',2)  
        plot([roof_x(4),roof_x(5)],[roof_y(4),roof_y(5)],'k','Linewidth',2)  
        plot([roof_x(5),roof_x(1)],[roof_y(5),roof_y(1)],'k','Linewidth',2)  
    else
        for i=1:Len
           Nod=El_Nod_Conn_Cover(i,:);
           plot([Global_Nodal_Coor_Cover(Nod(1),1),Global_Nodal_Coor_Cover(Nod(2),1)],[Global_Nodal_Coor_Cover(Nod(1),2),Global_Nodal_Coor_Cover(Nod(2),2)],'k','Linewidth',2)
           plot([Global_Nodal_Coor_Cover(Nod(2),1),Global_Nodal_Coor_Cover(Nod(3),1)],[Global_Nodal_Coor_Cover(Nod(2),2),Global_Nodal_Coor_Cover(Nod(3),2)],'k','Linewidth',2)   
           plot([Global_Nodal_Coor_Cover(Nod(3),1),Global_Nodal_Coor_Cover(Nod(4),1)],[Global_Nodal_Coor_Cover(Nod(3),2),Global_Nodal_Coor_Cover(Nod(4),2)],'k','Linewidth',2)
           plot([Global_Nodal_Coor_Cover(Nod(4),1),Global_Nodal_Coor_Cover(Nod(1),1)],[Global_Nodal_Coor_Cover(Nod(4),2),Global_Nodal_Coor_Cover(Nod(1),2)],'k','Linewidth',2)   
        end
    end
    xlim([-(Receiver_Width/2), (Receiver_Width/2)])
    ylim([0 Receiver_Width])
    title('Roof & Floor Discretization','FontSize', 28)
    xlabel('X position [m]','FontSize', 24) 
    ylabel('Y Position [m]','FontSize', 24) 
end

display('-----------------------------------------------')
display('Cavity Mesh Generation Completed')
display(['Panel Elements:' num2str(El_Panel) ' '])
display(['Roof/Floor Elements:' num2str(El_Cover) ' '])
display(['Total Number of Elements:' num2str(N_el_domain) ' '])
display('-----------------------------------------------')