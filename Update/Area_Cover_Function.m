function [A_Roof] = Area_Cover_Function(Receiver_Width,El_Nod_Conn_Cover,Global_Nodal_Coor_Cover)
    [Len,Wid]=size(El_Nod_Conn_Cover);
    for ii=1:Len %Element Loop
        for kk=1:4
            Coord(kk,1)=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,kk),1)]+(Receiver_Width/2);
            Coord(kk,2)=[Global_Nodal_Coor_Cover(El_Nod_Conn_Cover(ii,kk),2)];
        end
        Line_Check=zeros(1,4);

        jj=1; 
        Test_Pt(1)=Coord(1,2);
        Test_Pt(2)=Coord(2,2);
        Check=Test_Pt(2)-Test_Pt(1);
        if Check < 0 || Check > 0
            Line_Check(jj)=1;
        end

        jj=2; 
        Test_Pt(1)=Coord(2,1);
        Test_Pt(2)=Coord(3,1);
        Check=Test_Pt(2)-Test_Pt(1);
        if Check < 0 || Check > 0
            Line_Check(jj)=1;
        end

        jj=3; 
        Test_Pt(1)=Coord(3,2);
        Test_Pt(2)=Coord(4,2);
        Check=Test_Pt(2)-Test_Pt(1);
        if Check < 0 || Check > 0
            Line_Check(jj)=1;
        end

        jj=4; 
        Test_Pt(1)=Coord(4,1);
        Test_Pt(2)=Coord(1,1);
        Check=Test_Pt(2)-Test_Pt(1);
        if Check < 0 || Check > 0
            Line_Check(jj)=1;
        end

        if nnz(Line_Check)==2
            if Line_Check(3)==1 && Line_Check(4)==1
                Sq=(max(Coord(:,1))-min(Coord(:,1)))*(max(Coord(:,2))-min(Coord(:,2)));
                Tri_1=(abs(Coord(3,1)-Coord(4,1))*abs(Coord(3,2)-Coord(4,2)))/2;
                Tri_2=(abs(Coord(4,1)-Coord(1,1))*abs(Coord(4,2)-Coord(1,2)))/2;
                Sm_Sq=abs(Coord(3,2)-Coord(4,2))*abs(Coord(4,1)-Coord(1,1));
                A_Roof(ii)=Sq-Tri_1-Tri_2-Sm_Sq;
            elseif Line_Check(2)==1 && Line_Check(3)==1
                Sq=(max(Coord(:,1))-min(Coord(:,1)))*(max(Coord(:,2))-min(Coord(:,2)));
                Tri_1=(abs(Coord(3,1)-Coord(2,1))*abs(Coord(3,2)-Coord(2,2)))/2;
                Tri_2=(abs(Coord(4,1)-Coord(3,1))*abs(Coord(4,2)-Coord(3,2)))/2;
                Sm_Sq=abs(Coord(3,2)-Coord(4,2))*abs(Coord(2,1)-Coord(3,1));
                A_Roof(ii)=Sq-Tri_1-Tri_2-Sm_Sq;
            end
        elseif nnz(Line_Check)==1
            if Line_Check(2)==1
                Sq=(max(Coord(:,1))-min(Coord(:,1)))*(max(Coord(:,2))-min(Coord(:,2)));
                Tri=(abs(Coord(2,1)-Coord(3,1))*abs(Coord(2,2)-Coord(3,2)))/2;
                A_Roof(ii)=Sq-Tri;
            elseif Line_Check(3)==1
                Sq=(max(Coord(:,1))-min(Coord(:,1)))*(max(Coord(:,2))-min(Coord(:,2)));
                Tri=(abs(Coord(3,1)-Coord(4,1))*abs(Coord(3,2)-Coord(4,2)))/2;
                A_Roof(ii)=Sq-Tri;
            elseif Line_Check(4)==1
                Sq=(max(Coord(:,1))-min(Coord(:,1)))*(max(Coord(:,2))-min(Coord(:,2)));
                Tri=(abs(Coord(1,1)-Coord(4,1))*abs(Coord(1,2)-Coord(4,2)))/2;
                A_Roof(ii)=Sq-Tri;
            end
        else
            A_Roof(ii)=(max(Coord(:,1))-min(Coord(:,1)))*(max(Coord(:,2))-min(Coord(:,2)));
        end
    end
end

