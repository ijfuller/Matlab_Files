function [F_Panel_1_Roof,F_Panel_2_Roof] = VF_Reorder_Panel_Roof(panel_n,panel_m,cover_n,cover_m,F_Panel_1_Floor,F_Panel_2_Floor)
%Reorders matrix of panel-floor to panel-roof
    if panel_n==1
        F_Panel_1_Roof=flipud(F_Panel_1_Floor);
        F_Panel_2_Roof=flipud(F_Panel_2_Floor);
    else
        for i=1:(panel_m*panel_n)
            for j=1:(cover_n*cover_m)
                if i+(panel_n*panel_m)/2 > (panel_n*panel_m)
                    F_Panel_1_Roof(i,j)=F_Panel_1_Floor(i-(panel_n*panel_m)/2, j);
                    F_Panel_2_Roof(i,j)=F_Panel_2_Floor(i-(panel_n*panel_m)/2, j);
                else
                    F_Panel_1_Roof(i,j)=F_Panel_1_Floor(i+(panel_n*panel_m)/2, j);
                    F_Panel_2_Roof(i,j)=F_Panel_2_Floor(i+(panel_n*panel_m)/2, j);
                end
            end
        end
    end
end

