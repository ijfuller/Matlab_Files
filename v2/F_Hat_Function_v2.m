function [F_hat_sol,F_hat_therm,epsilon_sol,epsilon_therm,rho_sol,rho_therm] = F_Hat_Function_v2(N_el_domain,VF_Matrix, epsilon_sol,epsilon_therm)
    KD=eye(N_el_domain);
    
    if N_el_domain >= 8
        epsilon_sol(7)=epsilon_sol(5);
        epsilon_therm(7)=epsilon_therm(5);
        epsilon_sol(8)=1;
        epsilon_therm(8)=1;
    end
    if N_el_domain == 9
        epsilon_sol(8)=epsilon_sol(5);
        epsilon_therm(8)=epsilon_therm(5);
        epsilon_sol(9)=1;
        epsilon_therm(9)=1;
    end
    
    rho_sol=1-epsilon_sol; %Solar reflectivity vector
    F_hat_sol=inv(KD-VF_Matrix.*rho_sol')*VF_Matrix; %Solar F-Hat Matrix
    
    rho_therm=1-epsilon_therm;%Vector of thermal reflectivity
    F_hat_therm=inv(KD-VF_Matrix.*rho_therm')*VF_Matrix; %Thermal F-Hat Matrix
end