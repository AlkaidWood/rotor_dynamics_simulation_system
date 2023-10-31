clc
clear
close all
bar=waitbar(0,'建模中...');

%% Input initial parameters

InitialParameter = inputEssentialParameter();
InitialParameter = inputIntermediateBearing(InitialParameter );
InitialParameter = inputLoosingBearing(InitialParameter);
% InitialParameter = inputRubImpact(InitialParameter);

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


%% Input bifurcation parameters
bifurParNum = 500; 
bifurParStart = 0; %
bifurParEnd = 1000;
bifurPar = linspace(bifurParStart, bifurParEnd, bifurParNum);

%% Input calculate parameters
accelerationTime = 2; % s
stableTime = 3; % s
transientPeriodNum = 20;
findTime = 2; % find the maximum value in the last 2 seconds
vmax = 800;
acceleration = vmax / accelerationTime;
transientTime = transientPeriodNum*2*pi ./ vmax; % s
if vmax(1) == 0 % for rotating speed
    transientTime(1) = 0;
end

Parameter.Status.vmax = vmax;
Parameter.Status.acceleration = acceleration;

%% calcualte
% creat matrix save the result;
dofNum = Parameter.Mesh.dofNum;
bifurResultDis = cell(bifurParNum, 1);
bifurResultSpeed = cell(bifurParNum, 1);

% waite bar
str = '开始计算...';
waitbar(0,bar,str)
len = bifurParNum;

for iPar = 1:1:bifurParNum
    % change bifurcation parameter
    InitialParameter.LoosingBearing.loosingDamping = bifurPar(iPar);
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
    
    % find the index in t to match tSpan
    timeStart   = TEND - findTime; 
    timeEnd     = TEND;
    FINDERROR   = 0.00005;
    tStartIndex = find(( (timeStart-FINDERROR)<t & t<(timeStart+FINDERROR) ),1); 
    tEndIndex   = find(( (timeEnd-FINDERROR)<t & t<(timeEnd+FINDERROR) ),1); 
    
    %Initialize to save data
    saveDisplacement = cell(dofNum,1); 
    saveSpeed = cell(dofNum,1); 
    
    % calculate poincare point
    domega = Parameter.Status.vmax * [1; abs(Parameter.Status.ratio)]';
    Node = Parameter.Mesh.Node;
    dofOnNodeNo = Parameter.Mesh.dofOnNodeNo;
    dofOnShaftNo = [Node(dofOnNodeNo).onShaftNo];
    speed = 0;
    for iDof = 1:1:dofNum
        % check the speed
        if speed ~= domega(dofOnShaftNo(iDof))
            speed = domega(dofOnShaftNo(iDof));
            tCutPeriod = (2*pi)/speed;
            tCutPeriodIndex = floor(tCutPeriod*SAMPLINGFREQUENCY/REDUCEINTERVAL);
            tCutLoopNum = floor( (timeEnd-timeStart)/tCutPeriod )-1;
            tCutStartpoint = timeStart; % the start time for cutting signal (s) 
            tCutStartpointIndex = find(( (tCutStartpoint)<=t & t<(tCutStartpoint+0.0005) ),1); % find the index of time near the cut point
            [~, temporary] = max(  q(:,tCutStartpointIndex:(tCutStartpointIndex+tCutPeriodIndex)),[], 2 ); % find the maximum value near the start point
            temporary = temporary + tCutStartpointIndex - 1;
            tCutStartpointIndex = temporary; % the index of cuting time for each row
            tCutStartpoint = t(tCutStartpointIndex)'; % cut time for each row(displacement)
        end
        saveDisplacement{iDof} = zeros(tCutLoopNum,1); 
        saveSpeed{iDof} = zeros(tCutLoopNum,1); 
        for iLoop = 1:1:tCutLoopNum
            indexThisLoop = find((tCutStartpoint(iDof)+(iLoop-1)*tCutPeriod)<=t, 1);
            saveDisplacement{iDof}(iLoop) =  q(iDof,indexThisLoop);%t_cut_loop*dof
            saveSpeed{iDof}(iLoop) =  dq(iDof,indexThisLoop);%t_cut_loop*dof
        end % end for iLoop
    end % end for iDof
    
    % save poincare points
    bifurResultDis{iPar} = saveDisplacement;
    bifurResultSpeed{iPar} = saveSpeed;
    
    
    % waite bar
    str=['计算中...',num2str(100*iPar/len),'%'];
    waitbar(iPar/len,bar,str)
end

% save the result
save('bifurcationRotate', 'bifurPar', 'bifurResultDis', 'bifurResultSpeed', 'bifurParNum')

% close wait bar
delete(bar)


%% Plot
% figure name
dofNum = Parameter.Mesh.dofNum;
dofOnNodeNo = Parameter.Mesh.dofOnNodeNo;
figureIdentity  = cell(1,dofNum);
nodeNo          = 1;
dofInThisNode   = 0;
for iDof = 1:1:dofNum
    if nodeNo == dofOnNodeNo(iDof)
        dofInThisNode = dofInThisNode + 1;
    else
        nodeNo = nodeNo + 1;
        dofInThisNode = 1;
    end
    figureIdentity{iDof}=['Node-',num2str(dofOnNodeNo(iDof)),'-DOF-',num2str(dofInThisNode)];
end

% plot
dofNo = 1;
sz = 1; % pt
for iDof = 1:1:dofNum
     h = figure('name', figureIdentity{iDof});
    for iPar = 1:1:bifurParNum
        xspan = bifurPar(iPar)*ones(length(bifurResultDis{iPar}{iDof}),1);
        yspan = bifurResultDis{iPar}{iDof};
        scatter(xspan, yspan ,sz,...
                'filled','MarkerEdgeColor', [0.40784,0.5804,0.651],...
                'MarkerEdgeAlpha', .1,...
                'MarkerFaceColor', [0.40784,0.5804,0.651],...
                'MarkerFaceAlpha', .1); hold on
    end % end for iPar
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
        'FontSize'    , 7                           , ... 
        'FontName'    , 'Times New Roman'           ,...
        'xtick'       , linspace(0,1000,6)) 
    xlabel('$c_l$ (N$\cdot$s/m)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    isBearing = iDof>= 65 ;
    isTranslation = rem(iDof, 4)==1 || rem(iDof, 4)==2;
    if isBearing
        ylabel('$q$ (m)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    elseif isTranslation
        ylabel('$q$ (m)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    else
        ylabel('$q$ (rad)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
    end
    set(gcf,'Units','centimeters','Position',[6 6 7.2 5])
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    figureName =  figureIdentity{iDof};
    figurePath = 'G:/大学硕士/毕业论文/论文/result/loose/bifurcationDamping/';
    savefig(h,[figurePath, figureName, '.fig']);
    % print(h, [figurePath, figureName], '-depsc2')
    print(h, [figurePath, figureName], '-dpng','-r400')
    close(h)
end % end for iDof













