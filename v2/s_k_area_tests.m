clear all
clc

Receiver_Height=12;
Receiver_Width=14;

% Lip=[0:0.1:(Receiver_Height/2)];
Lip=2;
inclination=30;

for i=1:length(Lip)
    Top_Lip=Lip(i);
    Bottom_Lip=Lip(i);

    interior_angle=135; %Angle between two adjacent panels
    panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel

    Top_Lip_Area=Top_Lip*Receiver_Width;
    Bottom_Lip_Area=Bottom_Lip*Receiver_Width;
    Area_Roof=(1+sqrt(2))*panel_width^2;
    
    
    Inclination_Max_Panel=atand((Receiver_Height-Top_Lip)/(Receiver_Width/2))
    
    Inclination_Height=tand(inclination)*(Receiver_Width/2);
    
    
    A1=(panel_width*Receiver_Height*4)+(Area_Roof*2)+(Top_Lip_Area)+(Bottom_Lip_Area);
    A2=A1-(Bottom_Lip_Area);
    A3=A1-(Top_Lip_Area)-(Area_Roof)-(panel_width*Top_Lip*4);

    Mod(i)=(A1/A2)*(A3/A1)^0.63;
end
figure(1)
plot(Lip,Mod,'Linewidth',2)