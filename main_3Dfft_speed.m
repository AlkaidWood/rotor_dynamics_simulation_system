clc
clear
close all
bar=waitbar(0,'建模中...');

%% Input initial parameters

InitialParameter = inputEssentialParameter();
InitialParameter = inputIntermediateBearing(InitialParameter );
% InitialParameter = inputCouplingMisalignment(InitialParameter);
InitialParameter = inputLoosingBearing(InitialParameter);
% InitialParameter = inputRubImpact(InitialParameter);

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

cPar = [100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400];
cParNum = length(cPar);

%% Input calculate parameters
accelerationTime = 2; % s
stableTime = 3; % s
transientPeriodNum = 20;
findTime = 2; % find the maximum value in the last 2 seconds
vmax = cPar;
acceleration = vmax / accelerationTime;
transientTime = transientPeriodNum*2*pi ./ vmax; % s
if vmax(1) == 0 % for rotating speed
    transientTime(1) = 0;
end


%% Calculate
% creat matrix save the result;
dofNum = Parameter.Mesh.dofNum;
resultDis = cell(cParNum, 1);


% waite bar
str = '开始计算...';
waitbar(0,bar,str)
len = cParNum;

for iPar = 1:1:cParNum
    % change running statues
    Parameter.Status.vmax = vmax(iPar);
    Parameter.Status.acceleration = acceleration(iPar);
                           
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
    fL = vmax; % rad
    fH = vmax*abs(Parameter.Status.ratio);
    dataNum = 20;
    x1 = fL;
    x2 = fH;
    
    plot(x1/2/pi,x1, 'k--', x2/2/pi,x1, 'k--'); 
    
    

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
    yMax = max(cPar);
    yMin = min(cPar);
    xMax = 600;
    xlim([0 xMax])
    ylim([yMin-(yMax-yMin)*0.1 yMax+(yMax-yMin)*0.1])
    zlabelname = '$|$P1$|$';
    xlabelname = '$f$ (Hz)';
    ylabelname = '$\dot{\Phi}$ (rad/s)';
    xlabel(xlabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    ylabel(ylabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    zlabel(zlabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    
    
    set(gcf,'Units','centimeters','Position',[6 6 7.2 5])
    
    % save figure
    figureName =  figureIdentity{idof};
    figurePath = 'G:/大学硕士/毕业论文/论文/result/loose/fftSpeed/';
    savefig(h,[figurePath, figureName, '.fig']);
    print(h, [figurePath, figureName], '-dpng','-r400')
    close(h)
    
end





