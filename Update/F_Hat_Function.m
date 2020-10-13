function [F_hat_sol,F_hat_therm,rho_sol,rho_therm] = F_Hat_Function(N_el_domain,VF_Matrix, epsilon_sol,epsilon_therm)
    KD=eye(N_el_domain);
    
    epsilon_sol=epsilon_sol';
    epsilon_therm=epsilon_therm';
    
    rho_sol=1-epsilon_sol; %Solar reflectivity vector
    F_hat_sol=inv(KD-VF_Matrix.*rho_sol')*VF_Matrix; %Solar F-Hat Matrix
    
    rho_therm=1-epsilon_therm;%Vector of thermal reflectivity
    F_hat_therm=inv(KD-VF_Matrix.*rho_therm')*VF_Matrix; %Thermal F-Hat Matrix
    
    display('-----------------------------------------------')
    display('F-hat Matrix Completed')
    display('-----------------------------------------------')
end