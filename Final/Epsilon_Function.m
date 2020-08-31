function [epsilon_sol,epsilon_therm] = Epsilon_Function(panel_n,panel_m,cover_n,cover_m,e_act_sol,e_pass_sol,e_act_therm,e_pass_therm,Top_Lip,Bottom_Lip)
    %%%Vector of solar absorbtivities
    for i=1:(4*panel_n*panel_m)
       epsilon_sol(i)=e_act_sol; 
    end
    for i=1+(4*panel_n*panel_m):(2*cover_n*cover_m)+(4*panel_n*panel_m)
        epsilon_sol(i)=e_pass_sol;
    end
    epsilon_sol(1+(2*cover_n*cover_m)+(4*panel_n*panel_m))=1;

    %%%Vector of thermal emissivities
    for i=1:(4*panel_n*panel_m)
       epsilon_therm(i)=e_act_therm; 
    end
    for i=1+(4*panel_n*panel_m):(2*cover_n*cover_m)+(4*panel_n*panel_m)
        epsilon_therm(i)=e_pass_therm;
    end
    epsilon_therm(1+(2*cover_n*cover_m)+(4*panel_n*panel_m))=1;

    if Top_Lip > 0
       epsilon_sol(1+(2*cover_n*cover_m)+(4*panel_n*panel_m))=e_pass_sol;
       epsilon_sol(2+(2*cover_n*cover_m)+(4*panel_n*panel_m))=1;
       epsilon_therm(1+(2*cover_n*cover_m)+(4*panel_n*panel_m))=e_pass_therm;
       epsilon_therm(2+(2*cover_n*cover_m)+(4*panel_n*panel_m))=1;
    end

    if Bottom_Lip > 0
       epsilon_sol(2+(2*cover_n*cover_m)+(4*panel_n*panel_m))=e_pass_sol;
       epsilon_sol(3+(2*cover_n*cover_m)+(4*panel_n*panel_m))=1;
       epsilon_therm(2+(2*cover_n*cover_m)+(4*panel_n*panel_m))=e_pass_therm;
       epsilon_therm(3+(2*cover_n*cover_m)+(4*panel_n*panel_m))=1;
    end
end

