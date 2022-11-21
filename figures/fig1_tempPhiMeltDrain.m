%% set up plot commands
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextFontSize',11)
addpath ../model_and_dependencies/
clear all; close all;
fs = 11;
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
fn2 = '03800'; % no drainge, Ea = 22 still ~0.5% drainage
inds2 = [2,2000,10000,17000,26000];

fn1 = '04304'; % partial drainage ~40%
inds1 = [2,1000,6600,9000,12500];

fns = {fn1, fn2};
inds = [inds1; inds2];

d = 10;

f = figure;
f.Units = 'centimeters';
% [left bottom width height]
f.Position = [1,1,21.5,9.5*2+0.5]; % should be 19 wide
axs = cell(2,10);

lMarg = 0.9;
bMarg = 1;
wid = 4;
hMarg = -0.2;
vMarg2 = 0.5;
vMarg1 = 0.25;
tFrame = zeros(2,5);
for i = 1:2
    for j = 1:5
        load([fp fns{i} '/i' num2str(inds(i,j)) '.mat'],'T','phi','tVec');
        TPlot = reshape(T,Grid.p.Ny,Grid.p.Nx)*173+100;
        phiPlot = reshape(phi,Grid.p.Ny,Grid.p.Nx);

%         TPlot(phiPlot > 1e-2) = nan;
        phiPlot(phiPlot < 1e-16) = nan;
        
        axs{i,j} = subplot(4,5,j+(i-1)*10);
        hold on
        contourf(X*d,Y*d,TPlot,40,'linestyle','none');
        caxis([100 273]);
        axis equal
        xlim([0 10]);
        ylim([-1 10]);
        axs{i,j}.XTickLabel = [];
        axs{i,j}.YTickLabel = [];
        axs{i,j}.Units = 'centimeters';
        colormap(axs{i,j},'bone');

        if i == 1 && j == 5
            c = colorbar('EastOutside');
            c.Label.Interpreter = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.Label.FontSize = fs;
            c.Label.String = 'Temperature, K';
            if i == 1
                cs{1} = c;
            else
                cs{2} = c;
            end
         end

        axs{i,5+j} = subplot(4,5,(i-1)*10+5+j);
        axs{i,5+j}.Colormap = colormap(flipud(parula));
        contourf(X*d,Y*d,phiPlot,40,'linestyle','none');
        axis equal
        xlim([0 10]);
        ylim([-1 10]);
        axs{i,5+j}.XTickLabel = [];
        axs{i,5+j}.YTickLabel = [];
        axs{i,5+j}.Units = 'centimeters';
        caxis([0 1]);
        
        
        if i == 1 && j == 5
            c = colorbar('EastOutside');
            c.Label.Interpreter = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.Label.FontSize = fs;%8;
            c.Label.String = 'Melt fraction';
            if i == 1
                cs{3} = c;
            else
                cs{4} = c;
            end
        end
        tFrame(i,j) = round(tVec(end),-2);
        
    end
end


%% rearrange
for i = 1:2
    for j = 1:5
        if i == 1
            axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+3*(wid+vMarg1)+vMarg2 , wid, wid];
        elseif i == 2
            axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+1*(wid+vMarg1), wid, wid];
        end
        
        if i == 1
            axs{i,5+j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+2*(wid+vMarg1)+vMarg2 , wid, wid];
        elseif i == 2
            axs{i,5+j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg, wid, wid];
        end
    end
end
        
yTikLab = {'0','5','10'};
xTikLab = {'','5','10'};
tik = [0,5,10];


% l b w h
i = 2;
j = 6;
% axs{i,j}.Position = [lMarg,bMarg,wid,wid];
axs{i,j}.YTick = tik;
axs{i,j}.YTickLabel = yTikLab;
axs{i,j}.XTick = tik;
axs{i,j}.XTickLabel = yTikLab;
axs{i,j}.YLabel.String = 'z-dir, km';
axs{i,j}.XLabel.String = 'radius, km';

axs{i,j}.YLabel.FontSize = fs;
axs{i,j}.XLabel.FontSize = fs;

i = 2;
j = 1;
% axs{i,j}.Position = [lMarg,bMarg+wid+vMarg1,wid,wid];
axs{i,j}.YTick = tik;
axs{i,j}.YTickLabel = yTikLab;
axs{i,j}.YLabel.String = 'z-dir, km';

axs{i,j}.YLabel.FontSize = fs;


i = 1;
j = 1;
% axs{i,j}.Position = [lMarg,bMarg+2*(wid+vMarg1)+vMarg2,wid,wid];
axs{i,j}.YTick = tik;
axs{i,j}.YTickLabel = yTikLab;
axs{i,j}.YLabel.String = 'z-dir, km';

axs{i,j}.YLabel.FontSize = fs;

i = 1;
j = 6;
% axs{i,j}.Position = [lMarg,bMarg+3*(wid+vMarg1)+vMarg2,wid,wid];
axs{i,j}.YTick = tik;
axs{i,j}.YTickLabel = yTikLab;
axs{i,j}.YLabel.String = 'z-dir, km';

axs{i,j}.YLabel.FontSize = fs;

for i = 1:length(cs)
    cs{i}.Position(3) = 0.01;
    cs{i}.Position(1) = 0.935;
end
%%
figLab = 'abcdefghij';
for i = 1:5
    colormap(axs{1,i},'bone'); 
    colormap(axs{2,i},'bone'); 
    set([axs{2,i},axs{2,i+5}],'Position',axs{1,i+5}.Position)
    axs{1,i+5}.Visible = 'off';
    axs{2,i+5}.Visible = 'off';
    set([axs{1,i},axs{1,i+5}],'Position',axs{1,i}.Position)
    text(axs{2,i+5},9.9,-0.55,[num2str(tFrame(2,i)) ' yrs'],...
        'FontSize',fs,'HorizontalAlignment','right','Color',[1-1e-12,1,1])
    text(axs{1,i+5},9.9,-0.55,[num2str(tFrame(1,i)) ' yrs'],...
        'FontSize',fs,'HorizontalAlignment','right','Color',[1-1e-12,1,1])
    text(axs{1,i+5},9.7,9.4,figLab(i),'fontweight','bold','interpreter','tex',...
        'FontSize',fs,'HorizontalAlignment','right','Color',[1-1e-12,1,1])
    text(axs{2,i+5},9.8,9.4,figLab(i+5),'fontweight','bold','interpreter','tex',...
        'FontSize',fs,'HorizontalAlignment','right','Color',[1-1e-12,1,1])
end


i = 2;
j = 1;
axs{i,j}.XTick = tik;
axs{i,j}.XTickLabel = yTikLab;
axs{i,j}.XLabel.String = 'radius, km';
axs{i,j}.XLabel.FontSize = fs;


i = 2;
j = 2;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;
axs{i,j}.XTickLabel = xTikLab;
axs{i,j}.XLabel.String = 'radius, km';
axs{i,j}.XLabel.FontSize = fs;


i = 2;
j = 3;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;
axs{i,j}.XTickLabel = xTikLab;
axs{i,j}.XLabel.String = 'radius, km';
axs{i,j}.XLabel.FontSize = fs;

i = 2;
j = 4;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;
axs{i,j}.XTickLabel = xTikLab;
axs{i,j}.XLabel.String = 'radius, km';
axs{i,j}.XLabel.FontSize = fs;

i = 2;
j = 5;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;
axs{i,j}.XTickLabel = xTikLab;
axs{i,j}.XLabel.String = 'radius, km';
axs{i,j}.XLabel.FontSize = fs;

i = 1;
j = 2;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;

i = 1;
j = 3;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;

i = 1;
j = 4;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;

i = 1;
j = 5;
axs{i,j}.XTick = tik;
axs{i,j}.YTick = tik;

% f.PaperPosition = [0.0177,4.6440,4.4646,5.7121];