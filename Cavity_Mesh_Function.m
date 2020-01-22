function [El_Nod_Conn_Panel, El_Nod_Conn_Cover, Global_Nodal_Coord_Panel, Global_Nodal_Coor_Cover, N_el_domain] = Cavity_Mesh_Function(Receiver_Height,Receiver_Width,panel_n,panel_m,cover_n,cover_m)

%%%User Inputs%%%
% Receiver_Height=12;  %Receiver opening height
% Receiver_Width=14;   %Reciever opening width
% 
% panel_n=1; % Number of elements in x axis for panel
% panel_m=10; % Number of elements in y axis for panel
% 
% cover_n=4; %Number of elements in x axis for top and bottom of cavity. %Some issue, needs to be even
% cover_m=3; %Number of elements in y axis for top and bottom of cavity



%%%%%%%%%%%%%Code Section%%%%%%%%%%%%%%%
interior_angle=135; %Angle between two adjacent panels
panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel
El_Panel=panel_n*panel_m; %Number of elements for a single panel
El_Cover=cover_n*cover_m; %Number of elements for the roof/floor
N_el_domain=4*(panel_n*panel_m)+2*(cover_n*cover_m)+1; %Number of elements in domain


%%%For panels
panel_x=[0:(panel_width/panel_n):panel_width]; %Panel nodal locations local in x axis
panel_y=[0:(Receiver_Height/panel_m):Receiver_Height]; %Panel nodal locations local in y axis

%%%Global Coordinates Panels%%%
Panels=4;
Global_Nodal_Coord_Panel=zeros(panel_n*panel_m,3);
x_coord=[0:(panel_width/panel_n):panel_width];
y_coord=[0:(Receiver_Height/panel_m):Receiver_Height];

k=1;
for i=1:panel_m+1
   for j=1:panel_n+1
      Global_Nodal_Coord_Panel(k,:,1)=[x_coord(j),y_coord(i),0];
      k=k+1;
   end
end

N_nod_x=panel_n+1;
N_nod_y=panel_m+1;
Nod_per_el=4; %Nodes per quadrilateral element
N_el=panel_n*panel_m; %Number of elements in domain for quadralateral
El_Nod_Conn_Panel=zeros(N_el,Nod_per_el); %Initializing the node connectivity matrix
i=1; 
for j=1:panel_m
    for k=1:panel_n
        i;
        El_Nod_Conn_Panel(i,:)=[k+(j-1)*N_nod_x, k+1+(j-1)*N_nod_x,k+N_nod_x+(j-1)*N_nod_x, k+N_nod_x+1+(j-1)*N_nod_x]; %Node connectivity for every element
        i=i+1;
    end
end

figure(1)
hold on
%%%plots nodal points of panel
for i=1:length(panel_x)
    for j=1:length(panel_y)
        scatter(panel_x(i),panel_y(j),'k','filled')
    end
end
xlim([0 panel_width])
ylim([0 Receiver_Height])

%%%plot lines of panel
for i=1:length(panel_x)
    line_x=ones(1,length(panel_y))*panel_x(i);
    plot(line_x,panel_y,'k')
end
for i=1:length(panel_y)
    line_y=ones(1,length(panel_x))*panel_y(i);
    plot(panel_x,line_y,'k')
end
title('Panel Discretization','FontSize', 28)
xlabel('X position [m]','FontSize', 24) 
ylabel('Y Position [m]','FontSize', 24) 
%%%Panel End

%%%Cavity Top & Bottom Cover
%%%Get four roof external points and plot them
panel_angle=[-(180-(interior_angle/2)-90)*3 -(180-(interior_angle/2)-90) (180-(interior_angle/2)-90) (180-(interior_angle/2)-90)*3]; %Calulates the angle rotation from centerline for panels
panel_radial_dist=cosd(panel_angle(2))*(Receiver_Width/2); %Distance between the cavity origin and middle of panel
roof_x=[(-Receiver_Width/2) (-Receiver_Width/2)*cosd(45) 0 (Receiver_Width/2)*cosd(45) (Receiver_Width/2)];
roof_y=[0 (Receiver_Width/2)*cosd(45) (Receiver_Width/2) (Receiver_Width/2)*cosd(45) 0];
figure(2)
hold on
scatter(roof_x,roof_y,'k','filled')
plot(roof_x,roof_y,'k','LineWidth',2)
xlim([-Receiver_Width/2 Receiver_Width/2])
ylim([0 Receiver_Width])
%%%Get four roof external points and plot them end
title('Top and Bottom Discretization','FontSize', 28)
xlabel('X position [m]','FontSize', 24) 
ylabel('Y Position [m]','FontSize', 24) 
    
if cover_m==1 && cover_n==1
    Global_Nodal_Coor_Cover=[roof_x;roof_y]';
    El_Nod_Conn_Cover=[1 5 4 3 2];
else
    %%%Top and Bottom Mesh%%%
    x_div_loc=[-Receiver_Width/2:(Receiver_Width/(cover_n)):Receiver_Width/2]; %Location for x divisions

    %which panel does each division intersect
    for i=1:(length(x_div_loc)-2)
        if roof_x(1)<=x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(2)
            Loc(i)=1;
        end
        if roof_x(2)<x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(3)
            Loc(i)=2;
        end
        if roof_x(3)<x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(4)
            Loc(i)=3;
        end
        if roof_x(4)<x_div_loc(i+1) && x_div_loc(i+1)<=roof_x(5)
            Loc(i)=4;
        end
    end

    %%%Equations for panels
    m(1)=(roof_y(1)-roof_y(2))/(roof_x(1)-roof_x(2));
    y_int(1)=roof_y(1)-m(1)*roof_x(1);
    m(2)=(roof_y(2)-roof_y(3))/(roof_x(2)-roof_x(3));
    y_int(2)=roof_y(2)-m(2)*roof_x(2);
    m(3)=(roof_y(3)-roof_y(4))/(roof_x(3)-roof_x(4));
    y_int(3)=roof_y(3)-m(3)*roof_x(3);
    m(4)=(roof_y(4)-roof_y(5))/(roof_x(4)-roof_x(5));
    y_int(4)=roof_y(4)-m(4)*roof_x(4);

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
    % y_cover_coord(3,cover_n)=y_cover_coord(2,cover_n);%temp fix
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


    %%%Plot Roof Descritization
    for i=1:length(Loc)
        plot([x_div_loc(i+1),x_div_loc(i+1)],[0,height(i)],'k')
    end
    for i=1:2:cover_n*2
        for j=2:cover_m
            plot([x_cover_coord(i,j),x_cover_coord(i+1,j)],[y_cover_coord(i,j),y_cover_coord(i+1,j)],'k')
        end
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
end


display('------------------------------------------')
display('Cavity Mesh Generation Completed')
display(['Panel Elements:' num2str(El_Panel) ' '])
display(['Roof/Floor Elements:' num2str(El_Cover) ' '])
display(['Total Number of Elements:' num2str(N_el_domain) ' '])
display('------------------------------------------')