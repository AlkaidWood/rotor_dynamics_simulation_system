% This program is to calcualte mode shape for rotor system

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
criticalFre = 80;
%criticalFre = [80, 131, 155, 526, 710, 1020, 1039]';
%criticalFre = linspace(1054,1060,25);
% --------------------------------------
dofNum = Parameter.Mesh.dofNum;
eigNum = 2*dofNum;
speedRatio = Parameter.Status.ratio;
exciteFreNum = length(criticalFre);
eigResult = zeros(eigNum, exciteFreNum);
exciteFre = criticalFre;
eigVector = zeros(eigNum, exciteFreNum);

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
    A=[-M\(C+G), -M\K;...
        eye(dofNum),zeros(dofNum,dofNum)];
    A = full(A);
    
    % calculate eigen vectors
    [V, D] = eig(A);
    eigResult(:, iFre) = diag(D); % calculate the eigenvalue
    
    % find the vector corrosponding to the critical speed
    [~, index] = min(abs(eigResult(:, iFre)-exciteFre(iFre)));
    trans = V(:,index); % complex number 
    eigVector(:,iFre) = real(trans) + imag(trans); 
    
    % waite bar
    str=['计算中...',num2str(100*iFre/len),'%'];
    waitbar(iFre/len,bar,str)
end
eigResult = imag(eigResult); 
eigResult = sort(eigResult); 

%% Calculate the node position
nodeDistance = Parameter.Mesh.nodeDistance;
shaftNum = length(nodeDistance);

%calculate the offset of each shafts due to the intermediate bearing
Shaft = Parameter.Shaft;
offsetPosition = zeros(Shaft.amount,1);
if isfield(Parameter,'IntermediateBearing')
InterBearing = Parameter.IntermediateBearing; % short the variable
for iShaft = 1:1:Shaft.amount
if iShaft ==1
    offsetPosition(1) = 0; % position of shaft 1 is reference
else
    for iInterBearing = 1:1:InterBearing.amount
        if InterBearing.betweenShaftNo(iInterBearing,2) == iShaft
           basicShaftD = InterBearing.positionOnShaftDistance(iInterBearing, 1);
           laterShaftD = InterBearing.positionOnShaftDistance(iInterBearing, 2);
           basicShaftNo = InterBearing.betweenShaftNo(iInterBearing,1);
           offsetPosition(iShaft) = basicShaftD - laterShaftD...
                                    + offsetPosition(basicShaftNo);
        else
           offsetPosition(iShaft) = 0;
        end % if InterBearing.betweenShaftNo(iInterBearing,2) == iShaft
    end % for iInterBearing = 1:1:InterBearing.amount
end % if iShaft ==1
end % for iShaft = 1:1:Shaft.amount 
end % if isfiled(InitialParameter,'IntermediateBearing')

% offset shaft
for iShaft = 1:1:shaftNum
    nodeDistance{iShaft} = nodeDistance{iShaft} + offsetPosition(iShaft);
end

%% Extract displacement
direction = 1; % 1-x; 2-y; 3-\theta_x; 4-\theta_y
% get dof index in setting direction
targetDof = zeros(length(Node),1);
counter = 1;
for iNode=1:1:length(Node)
    if ~Node(iNode).isBearing
        targetDof(counter) = Parameter.Mesh.dofInterval(iNode,1)+direction-1;
        counter = counter + 1;
    end
end
targetDof(targetDof==0) = [];
eigVectorPart = eigVector(targetDof,:);

% divided eigen vector according to shaft No.
eigVectorShaft = cell(shaftNum,1);
counter = 1;
for iShaft = 1:1:shaftNum
    nodeNumThisShaft = length(nodeDistance{iShaft});
    eigVectorShaft{iShaft} = eigVectorPart(counter:counter+nodeNumThisShaft-1, :);
    counter = counter + nodeNumThisShaft;
end

%% Plot the mode shape
close all
mycolor = [0.94902,0.37647,0.32157;...
           0.40784,0.5804,0.651];
for iFre = 1:1:exciteFreNum
    h = figure;
    for iShaft = 1:1:shaftNum
        distance = nodeDistance{iShaft};
        shape = eigVectorShaft{iShaft}(:, iFre);
        c = polyfit(distance, shape, 3); %进行拟合，c为2次拟合后的系数
        meshValue = linspace(distance(1),distance(end),100);
        d = polyval(c, meshValue, 1); %拟合后，每一个横坐标对应的值即为d
        plot(meshValue, d,...
             'color', mycolor(iShaft, :),...
             'linewidth',1.5); hold on
        scatter(distance, shape, 24,...
        'filled',...
        'MarkerFaceColor',mycolor(iShaft,:),...
        'HandleVisibility','off');
    end
    hold off
    ylim([-0.4, 0.4])
    xlim([-0.2,2.6])
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
    legend('低压转子','高压转子', 'Fontname', '宋体', 'FontSize',8);
    xlabel('\fontname{宋体}节点位置 \fontname{Times New Roman}(m)','FontSize',9);
    ylabel('振型','Fontname', '宋体', 'FontSize',9);
    set(gcf,'Units','centimeters','Position',[6 6 14 4.5])
    figureName =  'fullModelShape';
    %figurePath = 'modeShape/';
    figurePath = 'G:/大学硕士/毕业论文/论文/result/fullModel/固有特性/';
    savefig(h,[figurePath, figureName, '.fig']);
    print(h, [figurePath, figureName], '-depsc2');
end





