%%%Monte_Carlo roof to floor
function [F_Roof_Floor] = MC_Roof_to_Floor_Pent_Function(Receiver_Height,Receiver_Width,N_rays,Global_Nodal_Coor_Cover)

plot_rays=0; %Disable=0, Enable=1.

x=Global_Nodal_Coor_Cover(:,1)+Receiver_Width/2;
y=Global_Nodal_Coor_Cover(:,2);
z=ones(1,5)*Receiver_Height;

if plot_rays==1
    figure(1)
    hold on
    grid on
    xlabel('X') 
    ylabel('Y')
    zlabel('Z')
    xlim([-14, 14*2])
    ylim([-14, 14*2])
    zlim([0, 12*2])

    scatter3(x,y,zeros(1,5),'filled','k')
    scatter3(x,y,z,'filled','b')

    plot3([x(1) x(2)],[y(1) y(2)],[0 0],'k','LineWidth',2)
    plot3([x(2) x(3)],[y(2) y(3)],[0 0],'k','LineWidth',2)
    plot3([x(3) x(4)],[y(3) y(4)],[0 0],'k','LineWidth',2)
    plot3([x(4) x(5)],[y(4) y(5)],[0 0],'k','LineWidth',2)
    plot3([x(5) x(1)],[y(5) y(1)],[0 0],'k','LineWidth',2)

    plot3([x(1) x(2)],[y(1) y(2)],[z(1) z(2)],'b','LineWidth',2)
    plot3([x(2) x(3)],[y(2) y(3)],[z(2) z(3)],'b','LineWidth',2)
    plot3([x(3) x(4)],[y(3) y(4)],[z(3) z(4)],'b','LineWidth',2)
    plot3([x(4) x(5)],[y(4) y(5)],[z(4) z(5)],'b','LineWidth',2)
    plot3([x(5) x(1)],[y(5) y(1)],[z(5) z(1)],'b','LineWidth',2)
end


%%%Generate an origin point for ray
for ii=1
hits=0;
for i=1:N_rays
    Pass_Fail_Set=[0 0 0 0 0];
    while nnz(Pass_Fail_Set)<5
        x_random=rand*Receiver_Width;
        y_random=rand*Receiver_Width/2;
        [Pass_Fail_Set] = PI_Pent_Roof_Floor_Function(x_random,y_random);
    end
    if plot_rays==1
        scatter(x_random,y_random,'g','filled')
    end
    Ptheta=rand;
    theta=asin(sqrt(Ptheta));
    Pphi=rand;
    phi=Pphi*2*pi;

    rx=cos(phi)*sin(theta); %Unit x component
    ry=sin(phi)*sin(theta); %Unit y component
    rz=cos(theta); %Unit z component
    scal=abs(rz)/Receiver_Height; %Scalar value for unit vector 
    rx=rx/scal; %Intersection x component
    ry=ry/scal; %Intersection y component
    rz=rz/scal; %Intersection z component
    if plot_rays==1
        plot3([x_random,x_random+rx],[y_random,y_random+ry],[0,rz],'--k')
    end
        
    x_random=x_random+rx;
    y_random=y_random+ry;
    Pass_Fail_Set=[0 0 0 0 0];
    [Pass_Fail_Set] = PI_Pent_Roof_Floor_Function(x_random,y_random);
        if nnz(Pass_Fail_Set)==5
            if plot_rays==1
                scatter3(x_random,y_random,rz,'g','filled')
            end
            hits=hits+1;
        else
            if plot_rays==1
                scatter3(x_random,y_random,rz,'r','filled')
            end
        end
end
F_Roof_Floor(ii)=hits/N_rays;
end

