clc
close all
clear 

%% Input initial parameters

InitialParameter = inputEssentialParameter();
InitialParameter = inputIntermediateBearing(InitialParameter);


% If you need change some parameters, please change the data in the struct:
% InitialParameter, then use establishModel( ) to get the different models

%% Establish models

% grid{1} = [1 4 1]; manualGrid{2} = [3]; 
grid = 'low';
Parameter = establishModel(InitialParameter,...
                           'gridfineness', grid,...
                           'isPlotModel',  false,...
                           'isPlotMesh',   false);
save('modelParameter','Parameter')  

%% Input calculate parameters
t = linspace(0,10,1000);

%% calculate acceleration speed phase
% load constants
Status   = Parameter.Status;
shaftNum = Parameter.Shaft.amount;
dofNum   = Parameter.Mesh.dofNum;
vmax            = Status.vmax;
duration        = Status.duration;
acceleration    = Status.acceleration;
isDeceleration  = Status.isDeceleration;
vmin            = Status.vmin;
ratio           = Status.ratio;
ratio           = [1; ratio]; % the first shaft is basic


% initial the status matrix
omega 	= zeros(dofNum,length(t));
domega  = omega;
ddomega = omega;


% calculate the status matrix
status = cell(shaftNum,1);
for iShaft = 1:1:shaftNum
    iVmax           = vmax * ratio(iShaft);
    iVmin           = vmin * ratio(iShaft);
    iAcceleration   = acceleration*ratio(iShaft);
    status{iShaft} = rotationalStatus(t, iVmax, duration, iAcceleration,...
                                      isDeceleration, iVmin);
end

% assemble the status matrix
Node = Parameter.Mesh.Node;
dofOnNodeNo = Parameter.Mesh.dofOnNodeNo;
for iDof = 1:1:dofNum
    isBearing = Node( dofOnNodeNo(iDof) ).isBearing;
    if ~isBearing
        shaftNo = Node( dofOnNodeNo(iDof) ).onShaftNo;
        omega(iDof,:)   = status{shaftNo}(1,:);
        domega(iDof,:)  = status{shaftNo}(2,:);
        ddomega(iDof,:) = status{shaftNo}(3,:);
    end
end

%% Plot 
% line name
lineName = cell(length(status),1);
for iShaft = 1:1:length(status)
    lineName{iShaft} = ['Shaft ', num2str(iShaft)];
end
lineName{1} = '低压转轴';
lineName{2} = '高压转轴';

% define color
color_lin=[ 0.40784,0.5804,0.651;
            0.94902,0.37647,0.32157;
            0 0.30078125 0.62890625;
            0.8125 0.3125 0.36328125;
            0.51953125 0.29296875 0.59375;
            0.4765625 0.7734375 0.671875];
blockcolor1=[254, 239, 237]/255;
blockcolor2=[240, 244, 246]/255;

% plot speed
figureName = 'Speed in running status';
h2 = figure('name',figureName);
%add squre block(position)
ymax = 2500;
ymin = -2500;
t1 = 2.5;
t2 = 7.5;
t3 = 10;
xmax = t3;
x_point = [0 t1 t1 t2 t2 t3];
y_point = ymax;
xline = [x_point(1),x_point(1),x_point(2),x_point(2);
         x_point(3),x_point(3),x_point(4),x_point(4);
         x_point(5),x_point(5),x_point(6),x_point(6)];
y1ine = [ymin,y_point,y_point,ymin];
% plot block
fill_1=fill(xline(1,:),y1ine,blockcolor1, 'HandleVisibility','off');hold on
fill_2=fill(xline(2,:),y1ine,blockcolor2, 'HandleVisibility','off');hold on
fill_3=fill(xline(3,:),y1ine,blockcolor1, 'HandleVisibility','off');hold on
set(fill_1,'edgecolor','none');
set(fill_2,'edgecolor','none');
set(fill_3,'edgecolor','none');
% plot line
for iShaft = 1:1:length(status)
    plot(t,status{iShaft}(2,:),'linewidth', 1.5, 'color', color_lin(iShaft,:)); hold on
end
% set size
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , [0,0,0,0]                   , ...
    'TickDir'     , 'in'                        , ...
    'XMinorTick'  , 'off'                       , ...
    'YMinorTick'  , 'off'                       , ...
    'TickLength'  , [.01 .01]                   , ...
    'LineWidth'   , 0.5                         , ...
    'XGrid'       , 'off'                        , ...
    'YGrid'       , 'off'                        , ...
    'FontSize'    , 7                           ,...
    'layer'       , 'top') 
% plot label
fzt = 8;
texthight=1700;
text_1='瞬态加速'; text_2='稳态'; text_3='瞬态减速';
text((x_point(1)+x_point(2))/2, texthight,text_1,'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',fzt);
text((x_point(3)+x_point(4))/2, texthight,text_2,'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',fzt);
text((x_point(5)+x_point(6))/2, texthight,text_3,'HorizontalAlignment','center', 'Fontname', '宋体','FontSize',fzt);

legend(lineName{:},...
        'FontSize',8,...
        'location', 'NorthOutside',...
        'box', 'off',...
        'orientation','horizontal')
    
set(gcf, 'unit', 'centimeters', 'position', [22 12 14 5])
ylim([ymin ymax])
xlabel('$t$ (s)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
ylabel('$\dot{\Phi}$ (rad/s)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');

%% save figure
figureName2 =  'evolution';
figurePath = 'G:/大学硕士/毕业论文/论文/result/fullModel/瞬态响应/';
savefig(h2,[figurePath, figureName2, '.fig']);
print(h2, [figurePath, figureName2], '-depsc2')

