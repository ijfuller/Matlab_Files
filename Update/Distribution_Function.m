function [Q_nodal_perc] = Distribution_Function(panel_n,panel_m,Q_nodal_perc,Dist_Pattern)
    if Dist_Pattern==1
        Q_nodal_perc=ones(1,(panel_n*panel_m*4))*(1/(panel_n*panel_m));
    elseif Dist_Pattern==2
        [i,j]=size(Q_nodal_perc);
        if j > (panel_n*panel_m) || j < (panel_n*panel_m)
            display("Error! Q_nodal_perc dimensions do not match number of nodes.")
        end
        if sum(Q_nodal_perc) > 1 || sum(Q_nodal_perc) < 1
            display("Error! Q_nodal_perc does not add up to 100%.")
        end
        temp=Q_nodal_perc;
        for i=2:4
            temp=[temp Q_nodal_perc]; 
        end
        Q_nodal_perc=temp;
    else 
        [i,j]=size(Q_nodal_perc);
        if j > (panel_n*panel_m*4) || j < (panel_n*panel_m*4)
            display("Error! Q_nodal_perc dimensions do not match number of nodes.")
        end 
        if sum(Q_nodal_perc) > 1 || sum(Q_nodal_perc) < 1
            display("Error! Q_nodal_perc does not add up to 100%.")
        end
    end
end

