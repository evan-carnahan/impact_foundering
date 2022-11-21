%% make manannan drainage simulation
addpath ../model_and_dependencies/
clear all; close all;
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
fs = 11;
set(groot,'defaultTextFontSize',fs)
%% set up grid
grRes = 100;
ocTh = grRes/5;
Gridp.xmin = 0; Gridp.xmax = 2; Gridp.Nx = grRes; 
Gridp.ymin = -ocTh/grRes; Gridp.ymax = 1; Gridp.Ny = grRes+ocTh;
Gridp.geom = 'cylindrical_rz';
Grid = build_stokes_grid_cyl(Gridp);
[X,Y] = meshgrid(Grid.p.xc,Grid.p.yc);

%% make plot dimensions
fp = './';
fn = '03321';
% inds = [1,42,118,295,412];
inds = [20,200,300,500,1400];
d = 10;

tPlt = zeros(1,length(inds));
f = figure;
f.Units = 'centimeters';
% [left bottom width height]
f.Position = [1,1,19,13.5];
% f.Position = [3,3,19,6.1];
axs = cell(3,5);

im = imread('manannanCrater.png');
axs{1,1} = subplot(3,5,[1 2]);
imshow(im);


for j = 1:5
    load([fp fn '/i' num2str(inds(j)) '.mat'],'T','phi','tVec');
    TPlot = reshape(T,Grid.p.Ny,Grid.p.Nx)*173+100;
    phiPlot = reshape(phi,Grid.p.Ny,Grid.p.Nx);
    
    axs{2,j} = subplot(3,5,j+5);
    contourf(X*d,Y*d,TPlot,40,'linestyle','none');
    hold on
    contour(X*d,Y*d,phiPlot,'r','LevelList',5e-2);
    caxis([100 273]);
    axis equal
    xlim([0 10]);
    ylim([-1 10]);
    axs{2,j}.FontSize = fs;
    axs{2,j}.XLabel.FontSize = fs;
    axs{2,j}.YLabel.FontSize = fs;
    axs{2,j}.XTickLabel = [];
    axs{2,j}.YTickLabel = [];
    axs{2,j}.Units = 'centimeters';
    text(9.9,-0.5,[num2str(round(tVec(end),-2)) ' yrs'],'FontSize',8,'HorizontalAlignment','right')

    if j == 5
        c1 = colorbar('EastOutside');
        c1.Label.Interpreter = 'latex';
        c1.TickLabelInterpreter = 'latex';
        c1.Label.FontSize = fs;
        c1.Label.String = 'Temperature, K';
    end

    axs{3,j} = subplot(3,5,10+j);
    axs{3,j}.Colormap = colormap(parula);
    contourf(X*d,Y*d,phiPlot,40,'linestyle','none');
    axis equal
    xlim([0 10]);
    ylim([-1 10]);
    axs{3,j}.FontSize = fs;
    axs{3,j}.XLabel.FontSize = fs;
    axs{3,j}.YLabel.FontSize = fs;
    axs{3,j}.FontSize = fs;
    axs{3,j}.XTickLabel = [];
    axs{3,j}.YTickLabel = [];
    axs{3,j}.Units = 'centimeters';
    caxis([0 1]);
    text(9.9,-0.5,[num2str(round(tVec(end),-2)) ' yrs'],'FontSize',8,'HorizontalAlignment','right')

    if j == 5
        c2 = colorbar('EastOutside');
        c2.Label.Interpreter = 'latex';
        c2.TickLabelInterpreter = 'latex';
        c2.Label.FontSize = fs;
        c2.Label.String = 'Melt fraction';
    end

end

set(f,'defaultAxesColorOrder',[0,0,0;0,0,0]);
load([fp fn '/i13258.mat']);
axs{1,2} = subplot(3,5,[3 4 5]);
hold on
plot(tVec,phiDrain1Vec/phiOrig*100,'b');
plot(tVec(inds),phiDrain1Vec(inds)/phiOrig*100,'ob');
yyaxis right
plot(tVec,phiDrain1Vec/1e9,'b');
ylim([0 phiOrig/1e9])
xlim([0 2000]);
yyaxis left

%%
% position
lMarg = 0.9;
bMarg = 1.05;
wid = 3.5;
hMarg = -0.2;
vMarg2 = 0.2;
vMarg1 = 0.08;

% tiks
yTikLab = {'0','5','10'};
xTikLab = {'','5','10'};
tik = [0,5,10];

for i = 2:3
    for j = 1:5
%         axs{i,j}.XTick = [];
        if i == 2
            axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+(wid+vMarg1)+vMarg2 , wid, wid];
            axs{i,j}.XTick = tik;
%             axs{i,j}.XTickLabel = xTikLab;
        elseif i == 3
            axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg, wid, wid];
            axs{i,j}.XLabel.String = 'radius, km';
            axs{i,j}.XTick = tik;
            axs{i,j}.XTickLabel = xTikLab;
        end
        if j == 1
            axs{i,j}.YLabel.String = 'z-dir, km';
            axs{i,j}.YTick = tik;
            axs{i,j}.YTickLabel = yTikLab;
            if i == 3
                axs{i,j}.XTick = tik;
                axs{i,j}.XTickLabel = yTikLab;
            end
        end
    end
end

c1.Position(1) = 0.921;
c2.Position(1) = 0.921;

%% adjust top row dimensions
axs{1,1}.Units = 'centimeters';
% axs{1,1}.Position(1) = 0;
axs{1,1}.Position(1) = lMarg;
axs{1,1}.Position(2) = axs{2,1}.Position(2) + axs{2,1}.Position(4) + 0.2;
axs{1,1}.Position(3) = 7;
axs{1,1}.Position(4) = 5;


axs{1,2}.Units = 'centimeters';
% position assignement


axs{1,2}.XLabel.String = 'time since impact, yrs';
axs{1,2}.FontSize = fs;
y1 = ylabel('Melt chamber drained, \%','FontSize',fs);
yyaxis right
y2 = ylabel('Melt volume drained, km$^3$','FontSize',fs);
axs{1,2}.Box = 'On';
axs{1,2}.Position = [9.5, 9.57, 8.139, 3.48];
%% rearrange - [l b w h]

%% add labels
figLab = 'cdefghijkl';
for i = 1:5
    text(axs{2,i},9.7,9.4,figLab(i),'fontweight','bold','interpreter','tex',...
        'FontSize',fs,'HorizontalAlignment','right','Color',[1-1e-12,1,1])
    text(axs{3,i},9.8,9.4,figLab(i+5),'fontweight','bold','interpreter','tex',...
        'FontSize',fs,'HorizontalAlignment','right','Color',[1-1e-12,1,1])
end

annotation('textbox',[0.06 0.97 0.1 0.1],...
'string','a','fontsize',fs,...
'fontweight','bold','interpreter','tex',...
'HorizontalAlignment','left','VerticalAlignment','bottom',...
'Margin',0,'LineStyle','none')

annotation('textbox',[0.415 0.97 0.1 0.1],...
'string','b','fontsize',fs,...
'fontweight','bold','interpreter','tex',...
'HorizontalAlignment','left','VerticalAlignment','bottom',...
'Margin',0,'LineStyle','none')
