clear all
clc

Receiver_Width=14;
interior_angle=135; %Angle between two adjacent panels
panel_width=(Receiver_Width/2)*cosd(interior_angle/2)*2; %Width of each panel


(((Receiver_Width/2)/sind(22.5))*sind(90))-panel_width
(((Receiver_Width/2)/sind(22.5))*sind(67.5))-(Receiver_Width/2)
(((Receiver_Width/2)/sind(22.5))*sind(67.5))