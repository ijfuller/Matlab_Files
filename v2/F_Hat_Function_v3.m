function [F_hat_sol,F_hat_therm,epsilon_sol,epsilon_therm,rho_sol,rho_therm] = F_Hat_Function_v3(N_el_domain,VF_Matrix, epsilon_sol,epsilon_therm,Top_Lip,Bottom_Lip,e_pass_sol,e_pass_therm,panel_m)
    KD=eye(N_el_domain);
    
    epsilon_sol=epsilon_sol';
    epsilon_therm=epsilon_therm';
    
    if Top_Lip > 0
       epsilon_sol(3+4*panel_m)=e_pass_sol;
       epsilon_sol(4+4*panel_m)=1;
       epsilon_therm(3+4*panel_m)=e_pass_therm;
       epsilon_therm(4+4*panel_m)=1;
    end
    

    if Bottom_Lip > 0
       epsilon_sol(4+4*panel_m)=e_pass_sol;
       epsilon_sol(5+4*panel_m)=1;
       epsilon_therm(4+4*panel_m)=e_pass_therm;
       epsilon_therm(5+4*panel_m)=1;
    end
    

    
    rho_sol=1-epsilon_sol; %Solar reflectivity vector
    F_hat_sol=inv(KD-VF_Matrix.*rho_sol')*VF_Matrix; %Solar F-Hat Matrix
    
    rho_therm=1-epsilon_therm;%Vector of thermal reflectivity
    F_hat_therm=inv(KD-VF_Matrix.*rho_therm')*VF_Matrix; %Thermal F-Hat Matrix
end