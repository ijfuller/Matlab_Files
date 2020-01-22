% function [outputArg1,outputArg2] = View_Factor_Matrix_Function(N_el_domain)
% %UNTITLED3 Summary of this function goes here
% %   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
% end

clear all
clc

T_amb=300;
T_Panel=350;
T_FR=300;

N_el_domain=7;
VF_Matrix=zeros(N_el_domain,N_el_domain); %Initialize fiew factor full matrix
rho=ones(7,1)*0.01; %Reflectivity of each element
epsilon=1-rho; %Emissivity of each element
Temp=ones(7,1)*300; %Temperature in K
for i=1:4
    Temp(i)=T_Panel;
end

A_1=5.38*12;%Panel_Area
A_2=69.297; %Floor Area
A_3=12*14; %Aperature Area

A(1)=A_1;
A(2)=A_1;
A(3)=A_1;
A(4)=A_1;
A(5)=A_2;
A(6)=A_2;
A(7)=A_3;

%Panel 1
VF_Matrix(1,1)=0; %Panel 1 to panel 1
VF_Matrix(1,2)=0.06138; %Panel 1 to panel 2
VF_Matrix(1,3)=0.08354; %Panel 1 to panel 3
VF_Matrix(1,4)=0.0918; %Panel 1 to panel 4
VF_Matrix(1,7)=0.4822; %Panel 1 to opening
VF_Matrix(1,5)=(1-sum(VF_Matrix(1,:)))/2; %Panel 1 to floor
VF_Matrix(1,6)=VF_Matrix(1,5); %Panel 1 to roof

%Panel 2
VF_Matrix(2,1)=0.06138; %Panel 2 to panel 1
VF_Matrix(2,2)=0; %Panel 2 to panel 2
VF_Matrix(2,3)=0.06138; %Panel 2 to panel 3
VF_Matrix(2,4)=0.08354; %Panel 2 to panel 4
VF_Matrix(2,7)=0.4717; %Panel 2 to opening
VF_Matrix(2,5)=(1-sum(VF_Matrix(2,:)))/2; %Panel 2 to floor
VF_Matrix(2,6)=VF_Matrix(2,5); %Panel 2 to roof

%Panel 3
VF_Matrix(3,1)=0.08354; %Panel 3 to panel 1
VF_Matrix(3,2)=0.06138; %Panel 3 to panel 2
VF_Matrix(3,3)=0; %Panel 3 to panel 3
VF_Matrix(3,4)=0.06138; %Panel 3 to panel 4
VF_Matrix(3,7)=0.4717; %Panel 3 to opening
VF_Matrix(3,5)=(1-sum(VF_Matrix(3,:)))/2; %Panel 3 to floor
VF_Matrix(3,6)=VF_Matrix(3,5); %Panel 3 to roof

%Panel 4
VF_Matrix(4,1)=0.0918; %Panel 4 to panel 1
VF_Matrix(4,2)=0.08354; %Panel 4 to panel 2
VF_Matrix(4,3)=0.06138; %Panel 4 to panel 3
VF_Matrix(4,4)=0; %Panel 4 to panel 4
VF_Matrix(4,7)=0.4822; %Panel 4 to opening
VF_Matrix(4,5)=(1-sum(VF_Matrix(4,:)))/2; %Panel 4 to floor
VF_Matrix(4,6)=VF_Matrix(4,5); %Panel 4 to roof

%Floor
VF_Matrix(5,1)=(VF_Matrix(1,5)*A_1)/A_2; %Floor to panel 1
VF_Matrix(5,2)=(VF_Matrix(2,5)*A_1)/A_2; %Floor to panel 2
VF_Matrix(5,3)=(VF_Matrix(3,5)*A_1)/A_2; %Floor to panel 3
VF_Matrix(5,4)=(VF_Matrix(4,5)*A_1)/A_2; %Floor to panel 4
VF_Matrix(5,7)=(0.1350*A_3)/A_2; %Floor to opening
VF_Matrix(5,5)=0; %Floor to floor
VF_Matrix(5,6)=1-sum(VF_Matrix(5,:)); %floor to roof

%Roof
VF_Matrix(6,:)=VF_Matrix(5,:);
VF_Matrix(6,5)=VF_Matrix(6,6);
VF_Matrix(6,6)=0;

%Aperature
VF_Matrix(7,1)=(VF_Matrix(1,7)*A_1)/A_3; %%Aperature to panel 1
VF_Matrix(7,2)=(VF_Matrix(2,7)*A_1)/A_3; %%Aperature to panel 2
VF_Matrix(7,3)=(VF_Matrix(3,7)*A_1)/A_3; %%Aperature to panel 3
VF_Matrix(7,4)=(VF_Matrix(4,7)*A_1)/A_3; %%Aperature to panel 4
VF_Matrix(7,7)=0; %%Aperature to opening
VF_Matrix(7,5)=(1-sum(VF_Matrix(7,:)))/2; %%Aperature to floor
VF_Matrix(7,6)=VF_Matrix(7,5); %%Aperature to roof

%E=(Temp.^4).*(5.67*10^-8);

for i=1:7
    for j=1:7
        [i,j]
        if i==j
            K_Delta=1
        else
            K_Delta=0
        end
        E(i,j)=(Temp(i)^4-Temp(j)^4)*(5.67*10^-8);
        F_hat(i,j)=inv((K_Delta-E(i,j)*rho(j)))*E(i,j)
    end
end