function [F_Roof_Floor] = MC_Roof_to_Floor_Function(Receiver_Width,Receiver_Height,cover_n,cover_m,N_rays,El_Nod_Conn_Cover,Global_Nodal_Coor_Cover)

    plot_MC_Ray=0;

    if plot_MC_Ray==1
        figure(4) 
        hold on
        grid on
        xlabel('X Axis') 
        ylabel('Y Axis') 
        zlabel('Z Axis')
        xlim([-(Receiver_Width/2),(Receiver_Width/2)]);
        ylim([0,Receiver_Width/2]);
        zlim([0,Receiver_Height]);
    end

    [Len,Wid]=size(El_Nod_Conn_Cover);
    for ii=1:Len
        for jj=1:Len
            %%%Coordinates of bottom panel
            x1=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,1),1);
            x2=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,2),1);
            x3=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,3),1);
            x4=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,4),1);
            y1=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,1),2);
            y2=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,2),2);
            y3=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,3),2);
            y4=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,4),2);
            %%%Coordinates of bottom panel
            x5=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),1);
            x6=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),1);
            x7=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),1);
            x8=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),1);
            y5=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,1),2);
            y6=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,2),2);
            y7=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,3),2);
            y8=Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(jj,4),2);

            if plot_MC_Ray==1
                plot([x1 x2],[y1 y2],'k','LineWidth',2)
                plot([x2 x3],[y2 y3],'k','LineWidth',2)
                plot([x3 x4],[y3 y4],'k','LineWidth',2)
                plot([x1 x4],[y1 y4],'k','LineWidth',2)

                plot3([x5 x6],[y5 y6],[Receiver_Height,Receiver_Height],'b','LineWidth',2)
                plot3([x6 x7],[y6 y7],[Receiver_Height,Receiver_Height],'b','LineWidth',2)
                plot3([x7 x8],[y7 y8],[Receiver_Height,Receiver_Height],'b','LineWidth',2)
                plot3([x8 x5],[y8 y5],[Receiver_Height,Receiver_Height],'b','LineWidth',2)
            end

            %%%Monte Carlo%%%
            hits=0;
            for i=1:N_rays
                x_cent=(x1+x2+x3+x4)/4;
                y_cent=(y1+y2+y3+y4)/4;

                m=[(y2-y1)/(x2-x1),(y3-y2)/(x3-x2),(y4-y3)/(x4-x3),(y1-y4)/(x1-x4)];
                y_int=[y1-m(1)*x1,y2-m(2)*x2,y3-m(3)*x3,y4-m(4)*x4];

                Pass_Fail_Set=zeros(4,1);
                x=[x1,x2,x3,x4];
                y=[y1,y2,y3,y4];
                while nnz(Pass_Fail_Set)<4
                    Pass_Fail_Set=zeros(4,1);
                    if cover_m==1
                        x_rand=rand*(Receiver_Width/(cover_n-0.25))-(Receiver_Width/(cover_n-0.25))/2+x_cent;
                        y_rand=rand*((Receiver_Width/2)/(cover_m-0.25))-((Receiver_Width/2)/(cover_m-0.25))/2+y_cent;
                    else
                        x_rand=rand*(Receiver_Width/(cover_n-1))-(Receiver_Width/(cover_n-1))/2+x_cent;
                        y_rand=rand*((Receiver_Width/2)/(cover_m-1))-((Receiver_Width/2)/(cover_m-1))/2+y_cent;
                    end
                    [Pass_Fail_Set] = PIIQ_0_Function(m,y_int,x_rand,y_rand,x,y);
                end
                if plot_MC_Ray==1
                    scatter(x_rand,y_rand,'filled','k')
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
                if plot_MC_Ray==1
                    plot3([x_rand,x_rand+rx],[y_rand,y_rand+ry],[0,rz],'--k')
                end

                Pass_Fail_Set=zeros(4,1);
                m=[(y6-y5)/(x6-x5),(y7-y6)/(x7-x6),(y8-y7)/(x8-x7),(y5-y8)/(x5-x8)];
                y_int=[y5-m(1)*x5,y6-m(2)*x6,y7-m(3)*x7,y8-m(4)*x8];
                x=[x5,x6,x7,x8];
                y=[y5,y6,y7,y8];
                x_rand=x_rand+rx;
                y_rand=y_rand+ry;
                [Pass_Fail_Set] = PIIQ_0_Function(m,y_int,x_rand,y_rand,x,y);
                if nnz(Pass_Fail_Set)==4
                    hits=hits+1;
                    if plot_MC_Ray==1 && nnz(Pass_Fail_Set)==4
                       scatter3(x_rand, y_rand,rz,'filled','g')
                    elseif plot_MC_Ray==1
                        scatter3(x_rand, y_rand,rz,'r')
                    end
                end
            end
            F(ii,jj)=hits/N_rays;
        end
    end
    F_Roof_Floor=F;

end

