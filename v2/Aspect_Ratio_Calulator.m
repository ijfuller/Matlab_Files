clear all
clc

display('---------------------')
%%%To find appropriate width
Receiver_Height=8
A_Panel=64.2908;
%%%%Calcs%%%
Panel_Width=A_Panel/Receiver_Height
interior_angle=135;
%panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel
Receiver_Width=(Panel_Width/(cosd(interior_angle/2)*2))*2
display('---------------------')

Aperature_Area=Receiver_Height*Receiver_Width
Roof_Area=(1+sqrt(2))*Panel_Width^2