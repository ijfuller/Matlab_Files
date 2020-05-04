%%%Works for 1D flux map
clear all
clc

% mu=0
% sigma=1
% 
% int_value=1;
% int=0.5;
% while int < 1
%     int_value=int_value+1;
%     x=[-int_value:0.1:int_value];
%     y=(exp(-(x-mu).^2/(2*sigma^2)))/(sigma*sqrt(2*pi));
% 
%     a=min(x);
%     b=max(x);
%     int=(erf((b-mu)/(sqrt(2)*sigma))-erf((a-mu)/(sqrt(2)*sigma)))/2;
% end

% y_nodes=10;
% x=[-int_value:((int_value*2)/y_nodes):int_value]
% for i=1:(length(x)-1)
%     a=x(i)
%     b=x(i+1)
%     Perc(i)=(erf((b-mu)/(sqrt(2)*sigma))-erf((a-mu)/(sqrt(2)*sigma)))/2
% end
% sum(Perc)


Q_in=50*10e6;
Q_in_per_panel=Q_in/4
Perc=[0.01 0.04 0.1 0.15 0.2];
Perc=[Perc fliplr(Perc)];


Receiver_Height=12;
y_nodes=10;
x_grid=[0,5.36];
y_grid=(-(Receiver_Height/2):(Receiver_Height/(y_nodes-1)):(Receiver_Height/2));

Perc=[Perc;Perc];

[X,Y]=meshgrid(x_grid,y_grid);
contourf(X,Y,Perc')
colorbar



