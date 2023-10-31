clc
clear
close all
bar=waitbar(0,'建模中...');

%% Input initial parameters

InitialParameter = inputEssentialParameter();
InitialParameter = inputIntermediateBearing(InitialParameter );
InitialParameter = inputCouplingMisalignment(InitialParameter);
% InitialParameter = inputLoosingBearing(InitialParameter);
InitialParameter = inputRubImpact(InitialParameter);

% If you need change some parameters, please change the data in the struct:
% InitialParameter, then use establishModel( ) to get the different models

%% Establish models

% grid{1} = [1 4 1]; manualGrid{2} = [3]; 
grid = 'low';
Parameter = establishModel(InitialParameter,...
                           'gridfineness', grid,...
                           'isPlotModel',  true,...
                           'isPlotMesh',   true);
save('modelParameter','Parameter')                       
%%  Generate the dynamic equation

generateDynamicEquation(Parameter);       

%% Input change parameter

cPar = linspace(0,7e6,8);
cParNum = length(cPar);

%% Input calculate parameters
accelerationTime = 2; % s
stableTime = 3; % s
transientPeriodNum = 20;
findTime = 2; % find the maximum value in the last 2 seconds
vmax = 263;
acceleration = vmax / accelerationTime;
transientTime = transientPeriodNum*2*pi ./ vmax; % s
if vmax(1) == 0 % for rotating speed
    transientTime(1) = 0;
end

InitialParameter.Status.vmax = vmax;
InitialParameter.Status.acceleration = acceleration;

%% Calculate
% creat matrix save the result;
dofNum = Parameter.Mesh.dofNum;
resultDis = cell(cParNum, 1);


% waite bar
str = '开始计算...';
waitbar(0,bar,str)
len = cParNum;

for iPar = 1:1:cParNum
    % change the parameter
    % InitialParameter.LoosingBearing.loosingStiffness = cPar(iPar);
    % InitialParameter.CouplingMisalignment.parallelOffset = cPar(iPar);
    InitialParameter.RubImpact.stiffness = cPar(iPar);
    grid = 'low';
    Parameter = establishModel(InitialParameter,...
                               'gridfineness', grid,...
                               'isPlotModel',  false,...
                               'isPlotMesh',   false);
                           
    % Generate the dynamic equation
    generateDynamicEquation(Parameter);
    
    % calculate response
    TSTART = 0;
    TEND = accelerationTime + transientTime + stableTime;
    SAMPLINGFREQUENCY = 10000;
    ISPLOTSTATUS = true;
    REDUCEINTERVAL = 1;
    tic
    [q, dq, t, convergenceStr] = calculateResponse(Parameter, [TSTART, TEND], SAMPLINGFREQUENCY,ISPLOTSTATUS,REDUCEINTERVAL);
    toc
    
    if ~isempty(convergenceStr)
        fprintf('%s \n', convergenceStr)
    end % end if 
    
    % save displacement
    resultDis{iPar} = q;

    % waite bar
    str=['计算中...',num2str(100*iPar/len),'%'];
    waitbar(iPar/len,bar,str)
end

%% Plot
tSpan = [2.5 5];
figureIdentity  = cell(1,dofNum);
nodeNo          = 1;
dofInThisNode   = 0;
dofOnNodeNo = Parameter.Mesh.dofOnNodeNo;
for idof = 1:1:dofNum
    if nodeNo == dofOnNodeNo(idof)
        dofInThisNode = dofInThisNode + 1;
    else
        nodeNo = nodeNo + 1;
        dofInThisNode = 1;
    end
    figureIdentity{idof}=['Node-',num2str(dofOnNodeNo(idof)),'-DOF-',num2str(dofInThisNode)];
end

% find the index in t to match tSpan
timeStart   = tSpan(1); 
timeEnd     = tSpan(2);
FINDERROR   = 0.00005;
tStartIndex = find(( (timeStart-FINDERROR)<t & t<(timeStart+FINDERROR) ),1); 
tEndIndex   = find(( (timeEnd-FINDERROR)<t & t<(timeEnd+FINDERROR) ),1); 

for idof = 1:1:dofNum
    
    h = figure('name', figureIdentity{idof});
    
    for iPar = 1:1:cParNum
        
        signal = resultDis{iPar}(idof,tStartIndex:tEndIndex);
        signallength = length(signal);

        % calculate
        Y  = fft(signal); 
        P2 = abs(Y/signallength); 
        P1 = P2(1:floor(signallength/2)+1);
        P1(2:end-1) = 2*P1(2:end-1); 
        f  = SAMPLINGFREQUENCY/REDUCEINTERVAL...
             *(0:(signallength/2))/signallength;
        xspan=f;
        yspan=P1;
        % plot
        plot3(xspan, cPar(iPar)*ones(length(xspan),1),yspan,'-','LineWidth',0.5,'color',[0 0.30078125 0.62890625]);hold on

    end
    % add lines
    fL = vmax/2/pi; % Hz
    fH = vmax*abs(Parameter.Status.ratio)/2/pi;
    dataNum = 20;
    y = linspace(cPar(1), cPar(end),dataNum);
    x1 = fL*ones(dataNum,1);
    x2 = 2 * fL*ones(dataNum,1);
    % x2 = fH*ones(dataNum,1);
    
    % plot(x1,y, 'k--'); 
    plot(x1,y, 'k--', x2,y, 'k--'); 
    
    % add text and arrow
    xMax = 300;
    xPos = xMax;
    axis = get(gca);
    yPos = abs((axis.YLim(2) - axis.YLim(1))) ;
    zPos = abs((axis.ZLim(2) - axis.ZLim(1)));
    text(xPos, axis.YLim(1)+ yPos*0.5, zPos*0.75,'碰摩刚度','rotation', -23,'VerticalAlignment','bottom', 'HorizontalAlignment','center','Fontname', '宋体','FontSize',8);
    text(xPos, axis.YLim(1)+ yPos*0.33, zPos*0.38,'轻微','rotation', -23,'VerticalAlignment','bottom', 'HorizontalAlignment','center','Fontname', '宋体','FontSize',8);
    text(xPos, axis.YLim(1)+ yPos*0.8, zPos*0.38,'严重','rotation', -23,'VerticalAlignment','bottom', 'HorizontalAlignment','center','Fontname', '宋体','FontSize',8);

    
    % 创建 arrow
annotation('arrow',[0.638765884652982 0.873765884652982],...
    [0.79450736546006 0.65666275393124],...
    'Color',[0.749019607843137 0.607843137254902 0.435294117647059],...
    'LineWidth',1.5,...
    'LineStyle',':',...
    'HeadWidth',5,...
    'HeadLength',5);

% 创建 doublearrow
annotation('doublearrow',[0.604166666666667 0.773958333333333],...
    [0.635382955771305 0.530744336569579],...
    'Color',[0.407843137254902 0.580392156862745 0.650980392156863],...
    'LineWidth',1.5,...
    'LineStyle',':',...
    'Head2Width',5,...
    'Head2Length',5,...
    'Head1Width',5,...
    'Head1Length',5);

% 创建 doublearrow
annotation('doublearrow',[0.783333333333333 0.9171875],...
    [0.522114347357066 0.442286947141316],...
    'Color',[0.949019607843137 0.376470588235294 0.32156862745098],...
    'LineWidth',1.5,...
    'LineStyle',':',...
    'Head2Width',5,...
    'Head2Length',5,...
    'Head1Width',5,...
    'Head1Length',5);

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
        'ZGrid'       , 'on'                        , ...
        'FontSize'    , 7                           , ... 
        'FontName'    , 'Times New Roman') 
    
    set(gca, 'YDir', 'reverse');
    
    xlim([0 xMax])
    zlabelname = '$|$P1$|$';
    xlabelname = '$f$ (Hz)';
    ylabelname = '$\Delta l$ (m)';
    xlabel(xlabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    ylabel(ylabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    zlabel(zlabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    
    
    set(gcf,'Units','centimeters','Position',[6 6 7.2 5])
    
    % save figure
    figureName =  figureIdentity{idof};
    figurePath = 'G:/大学硕士/毕业论文/论文/result/misalignmentRub/fftRubStiffness/';
    savefig(h,[figurePath, figureName, '.fig']);
    print(h, [figurePath, figureName], '-dpng','-r400')
    close(h)

end





