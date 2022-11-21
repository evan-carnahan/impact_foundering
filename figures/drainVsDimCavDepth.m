% make dimensionless penetration fraction vs. drainage plot
clear all; close all;
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextFontSize',11)
% impact keys
% icArr = {'03321','03800','04304','03314'...
%     ,'03402','03400','03313','03701','03330'};
Ht = [8.20,4.98,6.87,13,6.90,8.09,16.44,11.42,24];
d = [10,10,10,20,20,20,30,30,40];
drn = [0.93,0,0.41,0.72,0,0,0.72,0,0.62];

fs = 11;

dimPen = Ht./d;
lw = 2;
f = figure;
f.Units = 'centimeters';
f.Position = [1,1,9.5,11.5];
hold on;
set(gca,'FontSize',fs);
plot(dimPen(1:3),drn(1:3),'or','MarkerFaceColor','w','LineWidth',lw)
plot(dimPen(4:6),drn(4:6),'oc','MarkerFaceColor','w','LineWidth',lw)
plot(dimPen(7:8),drn(7:8),'ok','MarkerFaceColor','w','LineWidth',lw)
plot(dimPen(9),drn(9),'og','MarkerFaceColor','w','LineWidth',lw)
leg = legend('10 km','20 km','30 km','40 km','Location','NorthWest');
xlabel('Transient cavity depth/Ice shell thickness, $C/D$','FontSize',fs);
ylabel('Fraction of melt chamber drained','FontSize',fs);
xlim([0.3, 0.9]);
ax = f.CurrentAxes;
ax.Position(4) = ax.Position(3);
box on
% print(f,'/Users/evan/Documents/fig2_dimlessDrainage','-dpng','-r1000')