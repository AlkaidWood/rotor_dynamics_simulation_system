%% plot2DStandard
% plot a 2D line diagram with standerd setting
%% Syntax
% picture = plot2DStandard(x, y)
% 
% picture = plot2DStandard(x, y, xlabelname)
%
% picture = plot2DStandard(x, y, xlabelname, ylabelname)
%
% picture = plot2DStandard(x, y, xlabelname, ylabelname, isUsedInA4)
%
% picture = plot2DStandard(x, y, xlabelname, ylabelname, isUsedInA4, isOnlySet)
%% Description
% x, y: are vectors saving the data in x- y- direction
%
% xlabelname, ylabelname: are char 
%
% isUsedInA4: is a boolean controling the size of figure
%
% isOnlySet: is a boolean indicating weather plot( ) would be execute


function picture = plot2DStandard(x, y, xlabelname, ylabelname, isUsedInA4, isOnlySet)

% check input
if nargin < 6
    isOnlySet = false;
end

if nargin < 5
    isUsedInA4 = false;
end

if nargin < 4
    ylabelname = 'y';
end

if nargin < 3
    xlabelname = 'x';
end

%%

% plot
if ~isOnlySet
    picture = plot(x,y,'-','LineWidth',0.5,'color',[0 0.30078125 0.62890625]);%plot
else
    picture = [];
end

%%

% set 
legend off
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , [0,0,0,0]                   , ...
    'TickDir'     , 'in'                        , ...
    'XMinorTick'  , 'off'                       , ...
    'YMinorTick'  , 'off'                       , ...
    'TickLength'  , [.01 .01]                   , ...
    'LineWidth'   , 0.5                         , ...
    'XGrid'       , 'on'                        , ...
    'YGrid'       , 'on'                        , ...
    'FontSize'    , 7                          , ... 
    'FontName'    ,'Times New Roman'            ) 


if isUsedInA4
    set(gcf,'Units','centimeters','Position',[6 6 7.2 4]);%Set the size of figure(for A4)
else
    set(gcf,'Units','centimeters','Position',[6 6 10 5]);%Set the size of figure(for viewing)
end

xlabel(xlabelname, 'Interpreter','latex', 'Fontname', 'Times New Roman','FontSize',9);
ylabel(ylabelname, 'Interpreter','latex', 'Fontname', 'Times New Roman','FontSize',9);

end
