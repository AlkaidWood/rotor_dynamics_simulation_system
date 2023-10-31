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

%cPar = logspace(log10(2e6),log10(7e6),16);
cPar = linspace(6e-3,12e-3,16);
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
resultRubF = cell(cParNum, 1);

% paarmeter for rubimpact
rubDof = Parameter.Mesh.dofInterval(Parameter.RubImpact.positionOnShaftNode,:)';
rubDof = rubDof(1,:);
interval = Parameter.RubImpact.interval;
rubStiffness = Parameter.RubImpact.stiffness;
rubNum = length(rubDof);

% waite bar
str = '开始计算...';
waitbar(0,bar,str)
len = cParNum;

for iPar = 1:1:cParNum
    % change the parameter
    % InitialParameter.LoosingBearing.loosingStiffness = cPar(iPar);
    InitialParameter.CouplingMisalignment.parallelOffset = cPar(iPar);
    % InitialParameter.RubImpact.stiffness = cPar(iPar);
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
    
    % calculate rub force
    rubForce = zeros(rubNum, size(q,2));
    for iRub = 1:1:rubNum
        oneAmplitude = abs(sqrt(q(rubDof(iRub), :).^2 + q(rubDof(iRub)+1).^2));
        oneRubForce = zeros(1, size(q,2));
        for iTime = 1:1:size(q,2)
            if oneAmplitude(iTime) > interval
                oneRubForce(iTime) = (oneAmplitude(iTime) - interval(iRub)) * rubStiffness(iRub);
            end % end if
        end % end for iTime
        rubForce(iRub, :) = oneRubForce;
    end

    % save rub force
    resultRubF{iPar} = rubForce;

    % waite bar
    str=['计算中...',num2str(100*iPar/len),'%'];
    waitbar(iPar/len,bar,str)
end


%% Plot
tSpan = [2.5 5];
% find the index in t to match tSpan
timeStart   = tSpan(1); 
timeEnd     = tSpan(2);
FINDERROR   = 0.00005;
tStartIndex = find(( (timeStart-FINDERROR)<t & t<(timeStart+FINDERROR) ),1); 
tEndIndex   = find(( (timeEnd-FINDERROR)<t & t<(timeEnd+FINDERROR) ),1); 

rubNode = Parameter.RubImpact.positionOnShaftNode;
for iRub = 1:1:rubNum
    
    h = figure('name', ['Node-', num2str(rubNode(iRub))]);
    
    for iPar = 1:1:cParNum
        signal = resultRubF{iPar}(iRub,tStartIndex:tEndIndex);
        % signal(signal==0) = [];
        limRubForce = [max(signal); min(signal)];
        plot(cPar(iPar)*ones(2,1), limRubForce,...
             'LineStyle', ':',...
             'LineWidth', 1.5,...
             'Color', [0.749019607843137 0.607843137254902 0.435294117647059],...
             'Marker', 'square',...
             'MarkerSize', 5,...
             'MarkerEdgeColor', 'none',...
             'MarkerFaceColor', [0.749019607843137 0.607843137254902 0.435294117647059]); hold on
    end
    
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
    
    xlabelname = '$\Delta l$ (m)';
    ylabelname = '$F_n$ (N)';
    xlabel(xlabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    ylabel(ylabelname,'Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    
    set(gcf,'Units','centimeters','Position',[6 6 7.2 4])
    
    
    % save figure
    figureName =  ['Node-', num2str(rubNode(iRub))];
    figurePath = 'G:/大学硕士/毕业论文/论文/result/misalignmentRub/rubForceMisValue/';
    savefig(h,[figurePath, figureName, '.fig']);
    print(h, [figurePath, figureName, '.eps'], '-depsc2');
    print(h, [figurePath, figureName], '-dpng','-r400');
    close(h)
    
end





