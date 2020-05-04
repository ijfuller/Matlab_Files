%%%Works for 1D flux map
clear all
clc

mu=0
sigma=4

int_value=100;
% int=0.5;
% while int < 1
    x=[-int_value:0.1:int_value];
    y=(exp(-(x-mu).^2/(2*sigma^2)))/(sigma*sqrt(2*pi));

    a=min(x);
    b=max(x);
    int=(erf((b-mu)/(sqrt(2)*sigma))-erf((a-mu)/(sqrt(2)*sigma)))/2
    int_value=int_value+1;
% end

figure(1)
plot(x,y,'Linewidth',2)

Q_in=50*10e6;
A_Panel=64.2908;
Flux_per_panel=(Q_in/4)/A_Panel;
Receiver_Height=12;

nodes_panel_y=10;
dashes=[(-Receiver_Height/2):(Receiver_Height/nodes_panel_y):(Receiver_Height/2)];

x_mod=[(-Receiver_Height/2):((Receiver_Height/(length(x)-1))):(Receiver_Height/2)];
y_mod=y*Flux_per_panel;
figure(2)
hold on
plot(x_mod,y_mod,'Linewidth',2)
xlabel('Receiver Height Position') 
ylabel('Flux (MW/m^2)') 

max_value=max(y_mod);
for i=1:length(dashes)
    plot([dashes(i),dashes(i)],[0,max_value],'k--','Linewidth',1)
end

for i=1:(length(dashes)-1)
   a_mod=dashes(i);
   b_mod=dashes(i+1);
   Flux_Nodal(i)=Flux_per_panel*(erf((b_mod-mu)/(sqrt(2)*sigma))-erf((a_mod-mu)/(sqrt(2)*sigma)))/2;
   Power_Nodal=Flux_Nodal*(A_Panel/nodes_panel_y);
end
Total_Power=sum(Power_Nodal)

%%%Contour Plot of Power%%%
x_contour=[0,5.36];
y_contour=([-(Receiver_Height/2):(Receiver_Height/(length(Power_Nodal)-1)):(Receiver_Height/2)]);
[X,Y]=meshgrid(x_contour,y_contour);
Power_Nodal=[Power_Nodal; Power_Nodal];
figure(3)
hold on
contourf(X,Y,Power_Nodal')
colorbar
xlabel('X Position') 
ylabel('Y Position') 

for i=1:length(dashes)
    plot([0,5.36],[dashes(i),dashes(i)],'k--','Linewidth',1.5)
end
