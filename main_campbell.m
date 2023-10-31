% This program is to calcualte campbell for rotor system

clc
clear
close all
bar=waitbar(0,'建模中...');
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

%% Input calculating parameter
% --------------input-------------------
exciteFreStart = 0; % rad/s
exciteFreEnd = 1500; % rad/s
step = 10;
% --------------------------------------
dofNum = Parameter.Mesh.dofNum;
eigNum = 2*dofNum;
speedRatio = Parameter.Status.ratio;
exciteFreNum = (exciteFreEnd - exciteFreStart)/step + 1;
eigResult = zeros(eigNum, exciteFreNum);
exciteFre = linspace(exciteFreStart, exciteFreEnd, exciteFreNum);

%% Calculate eigenvalue 
% waite bar
str = '开始计算';
waitbar(0,bar,str)
len = exciteFreNum;

% get shaft dof range
shaftNum = Parameter.Shaft.amount;
Node = Parameter.Mesh.Node;
dofInterval = Parameter.Mesh.dofInterval;
shaftDof = zeros(shaftNum,2);
for iShaft = 1:1:shaftNum
    IShaftNode = Node( [Node.onShaftNo] == iShaft & [Node.isBearing] == false );
    startNode = min([IShaftNode.name]);
    endNode = max([IShaftNode.name]);
    shaftDof(iShaft,:) = [dofInterval(startNode,1), dofInterval(endNode,2)];
end

for iFre = 1:1:exciteFreNum
    % get matrix
    M = Parameter.Matrix.mass;
    C = Parameter.Matrix.damping;
    K = Parameter.Matrix.stiffness;
    G = Parameter.Matrix.gyroscopic;
   
    % update G with excite frequency
    iBasicSpeed = exciteFre(iFre);
    for iShaft = 1:1:shaftNum
        % get the spin speed for each shaft
        if iShaft==1
            shaftSpeed = iBasicSpeed;
        else
            shaftSpeed = shaftSpeed * speedRatio(iShaft-1);
        end
        G(shaftDof(iShaft,1):shaftDof(iShaft,2),shaftDof(iShaft,1):shaftDof(iShaft,2))...
            = shaftSpeed * G(shaftDof(iShaft,1):shaftDof(iShaft,2),shaftDof(iShaft,1):shaftDof(iShaft,2));
    end
    
    % assemble a new matrix A
    % A's eign is the eigenvalue of origin dydnamic system
    A=[-M\(C-G), -M\K;...
        eye(dofNum),zeros(dofNum,dofNum)];
    A = full(A);
    
    eigResult(:, iFre) = eig(A); % calculate the eigenvalue
    
    % waite bar
    str=['计算中...',num2str(100*iFre/len),'%'];
    waitbar(iFre/len,bar,str)
end
eigResult = imag(eigResult);
eigResult = sort(eigResult);

%% Filter
filterRatioMax = 2;
filterMaxValue = filterRatioMax * exciteFreEnd;
newEig = zeros(eigNum, exciteFreNum);
% delete lager value
counter = 1;
for iRow = 1:1:size(A,1)
    vector = abs(eigResult(iRow,:));
    isSmall = sum(find(vector<filterMaxValue));
    if isSmall
        newEig(counter,:) = eigResult(iRow,:);
        counter = counter + 1;
    end
end
newEig(all(newEig==0,2),:) = [];

% reshape the curve
nearPointsNum = 5;
trans = zeros(size(newEig,1,2));
for iColumn = 1:1:size(newEig,2)
    if iColumn==1
        trans(:, iColumn) = newEig(:, iColumn);
    elseif iColumn==2
        trans(:,iColumn) = newEig(:,iColumn);
        % calculate slope factor
        thisSlope = (trans(:, iColumn) - trans(:, iColumn-1)) / step;
    else
        % update data
        lastSlope = thisSlope;
        for iRow=1:1:size(newEig,1)
            % get the near points
            lastElement = trans(iRow, iColumn-1);
            distance = newEig(:,iColumn) - lastElement;
            [~, pointsIndex] = sort(abs(distance));
            points = newEig(pointsIndex,iColumn);
            points = points(1:nearPointsNum);
            % calculate slop factor
            pointsSlope = (points - lastElement) / step;
            [~, trueIndex] = min(abs(pointsSlope - lastSlope(iRow)));
            truePoint = points(trueIndex);
            truePoint = truePoint(1);
            trans(iRow,iColumn) = truePoint;
            thisSlope(iRow) = (trans(iRow, iColumn) - trans(iRow,iColumn-1)) / step;
        end
    end
end
newEig = trans;

% delete the curve changed under a ratio
changeRatio1 = 0.00; % (max - min)/exciteFreEnd < ratio
changeRatio2 = 5; % maxSlop>0 and minSlop<0 and ...
counter = 1;
trans = zeros(size(newEig,1,2));
for iRow = 1:1:size(newEig,1)
    thisSlop = zeros(size(newEig,2)-1,1);
    for iColumn = 2:1:size(newEig,2)
        thisSlop(iColumn-1) = (newEig(iRow,iColumn) - newEig(iRow,iColumn-1))/step;
    end
    condition1 = (max(newEig(iRow,:)) - min(newEig(iRow,:)))/exciteFreEnd <= changeRatio1;
    maxS = max(thisSlop);
    minS = min(thisSlop);
    condition21 = maxS >0 && minS<0;
    condition22 = (maxS - minS)/max(abs([maxS, minS]))> changeRatio2;
    condition2 = condition21 && condition22;
    isAbandon = condition1 || condition2;
    isSave = ~isAbandon;
    if isSave 
        trans(counter,:) = newEig(iRow,:);
        counter = counter + 1;
    end
end
newEig = trans;
newEig(all(newEig==0,2),:) = [];

%% Plot
h=figure;
for iCurve = 1:1:size(newEig,1)
    plot(exciteFre, newEig(iCurve, :),'linewidth',0.5,'color','black'); hold on
    %plot(exciteFre, trans(iCurve, :)); hold on
end
plot(exciteFre,exciteFre,'linewidth',1,'color',[0.94902,0.37647,0.32157]);
plot(exciteFre,1.3*exciteFre,'linewidth',1,'color',[0.40784,0.5804,0.651]);
%plot(exciteFre,-exciteFre,'linewidth',2);
%plot(exciteFre,-1.3*exciteFre,'linewidth',2);
node1= [80, 131, 155, 280, 375, 526, 558, 710, 868 1020, 1039, 1100, 1350];
node2= [60, 88, 122, 189, 300, 409, 420, 566, 660, 767, 824, 910, 1080, 1230, 1270];
scatter(node1, node1, 24,...
        'filled',...
        'MarkerFaceColor',[0.94902,0.37647,0.32157]); hold on
scatter(node2, 1.3*node2, 24,...
        'filled',...
        'MarkerFaceColor',[0.40784,0.5804,0.651]);
text(1290,1200, '$\Omega = \dot{\Phi}$',...
     'Fontname', 'Times New Roman',...
     'FontSize', 9,...
     'interpreter', 'latex')
 text(980,1200, '$\Omega = 1.3\dot{\Phi}$',...
     'Fontname', 'Times New Roman',...
     'FontSize', 9,...
     'interpreter', 'latex')
ylim([0 1750])

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
    'FontSize'    , 7                          , ... 
    'FontName'    ,'Times New Roman'            ) 
xlabel('$\dot{\Phi}$ (rad/s)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
ylabel('$\Omega$ (rad/s)','Fontname', 'Times New Roman', 'FontSize',9,'interpreter','latex');
set(gcf,'Units','centimeters','Position',[6 6 14 6])
figureName =  'fullModelCompbell';
%figurePath = 'G:/大学硕士/毕业论文/论文/result/fullModel/固有特性/';
%savefig(h,[figurePath, figureName, '.fig']);
%print(h, [figurePath, figureName], '-depsc2')



