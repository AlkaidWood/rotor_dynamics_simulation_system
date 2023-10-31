clc
close all
clear

%% Load the data
load postProcessData % get figureIdentity
targetPath = 'G:/大学硕士/毕业论文/论文/result/fullModel/分岔图/';
changePath= 'G:/大学硕士/毕业论文/论文/result/rub/分岔图/';

%% change the YLim of the bifurcation figure
figNum = length(figureIdentity);

for iFig = 1:1:figNum
    fullTargetPath = [targetPath, figureIdentity{iFig}];
    fullChangePath = [changePath, figureIdentity{iFig}];
    % get target YLim
    h = openfig(fullTargetPath);
    ha = get(gcf, 'Children');
    targetYLim = ha.YLim;
    close(h);
    
    % set the new YLim
    h = openfig(fullChangePath);
    ylim(targetYLim);
    savefig([fullChangePath,'.fig'])
    print(h, fullChangePath, '-dpng','-r400')
    close(h);
    
end



