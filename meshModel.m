%% inputEssentialParameter
% Meshing each shaft, output a struct including initial parameter and grid information
%% Syntax
% Parameter = meshModel(InitialParameter)
% Parameter = meshModel(InitialParameter,gridFineness) 
% Parameter = meshModel(InitialParameter,manualGrid)
%% Description
% Parameter: includes all information of InitialParameter and information
% about meshing result
%
% InitialParameter: generated by input... function
% 
% gridFineness = 'low' or 'middle' or 'high'
%
% manualGrid: is a 1*n cell (1*column). The grid information of the n-th
% shaft is saved in manualGrid{n}. 
% e.g. manualGrid{n} = [1 2 5 2] means there are four segments which are 
% divided into 1, 2, 5, 2 pieces，respectively.


function Parameter = meshModel(varargin)   
% define the constant: this function will merge key nodes whose distance
% from other nodes to itself < totalLength/THRESHOLD_COEFFICIENT
THRESHOLD_COEFFICIENT = 5000; 

% recognize the input parameter
if isempty(varargin)
    error('Insufficient input parameters')
else
    InitialParameter = varargin{1};
end

switch length(varargin)
    case 1
        isAutoMesh = true;
        gridFineness = 'low';
    case 2
        if ischar(varargin{2})
            gridFineness = varargin{2};
            isAutoMesh = true;
        elseif iscell(varargin{2})
            manualGrid = varargin{2};
            isAutoMesh =false;
        else
            error('Please input char or cell data for the second parameter')
        end
    otherwise 
        error('too much input parameter')
end

%%

% check the initial parameter
Shaft = InitialParameter.Shaft;
Disk = InitialParameter.Disk;
Bearing = InitialParameter.Bearing;

if InitialParameter.ComponentSwitch.hasRubImpact
    RubImpact = InitialParameter.RubImpact;
    hasRub = true;
else
    hasRub = false;
end

if InitialParameter.ComponentSwitch.hasIntermediateBearing
    InterBearing = InitialParameter.IntermediateBearing;
    hasInterBearing = true;
else
    hasInterBearing = false;
end

if InitialParameter.ComponentSwitch.hasCouplingMisalignment
    Coupling = InitialParameter.CouplingMisalignment;
    hasCoupling = true;
else
    hasCoupling = false;
end

if InitialParameter.ComponentSwitch.hasLoosingBearing
    LoosingBearing = InitialParameter.LoosingBearing;
    hasLoosingBearing = true;
else
    hasLoosingBearing = false;
end

%%

% record key points on the shaft
keyPoints = cell(Shaft.amount,1);

for iShaft = 1:1:Shaft.amount
    keyPoints{iShaft} = [0; Shaft.totalLength(iShaft)];
    
    
    % record disk
    position = find(Disk.inShaftNo == iShaft);
    keyPointsDisk = Disk.positionOnShaftDistance(position);
    
    
    % record bearing 
    position = find(Bearing.inShaftNo == iShaft);
    keyPointsBearing = Bearing.positionOnShaftDistance(position);
    
    % record rub
    if hasRub
        position = find(RubImpact.inShaftNo == iShaft);
        keyPointsRub = RubImpact.positionOnShaftDistance(position);
    else
        keyPointsRub = [];
    end
    
    
    % record intermediate bearing
    if hasInterBearing
        position = find(InterBearing.betweenShaftNo == iShaft);
        keyPointsInterBearing = InterBearing.positionOnShaftDistance(position);
    else
        keyPointsInterBearing = [];
    end
    
    
    % assembling
    keyPoints{iShaft} = sort([  keyPoints{iShaft};...
                                keyPointsDisk;...
                                keyPointsRub;...
                                keyPointsInterBearing;...
                                keyPointsBearing]);
          
                          
    % merge the same value (or too close) in keyPoints
    ii = 1;
    while ii<length(keyPoints{iShaft})
        distance = abs( keyPoints{iShaft}(ii) - keyPoints{iShaft}(ii+1) );
        isTooClose = distance < Shaft.totalLength(iShaft)/THRESHOLD_COEFFICIENT;
        if isTooClose
            keyPoints{iShaft}(ii+1) = [];
        else
            ii = ii + 1;
        end % end if
    end % end while

end % end for iShaft = 1:1:Shaft.amount

%%

% check the meshing parameter
if isAutoMesh
    
    switch gridFineness
        case 'low' 
            FINENESS = 1;
        case 'middle'
            FINENESS = 4;
        case 'high'
            FINENESS = 10;
        otherwise
            error('Please input: low, middle or high at the second parameter')
    end % end switch
    
else
    
    for iShaft = 1:1:Shaft.amount
        isMatch = length(manualGrid{iShaft}) == ( length(keyPoints{iShaft}) - 1 );
        if ~isMatch
            error(['the dimension of manual grid is not matched with',...
                   'segments. Please use auto mesh and check the number of segments'])
        end % end if
    end % end for iShaft
    
end % end if

%%

% mesh
rowSegmentNum = zeros(Shaft.amount,1);
nodeDistance = cell(Shaft.amount,1);
for iShaft = 1:1:Shaft.amount
    % calculate the number of segments for each shaft
    rowSegmentNum(iShaft) = length(keyPoints{iShaft}) - 1; 
    
    if isAutoMesh
        % calculate the number of elements for each segment
        standardLength = Shaft.totalLength(iShaft)/ FINENESS;
        elementNum = ceil( diff(keyPoints{iShaft}) ./ standardLength );
    else
        elementNum = manualGrid{iShaft};
    end
    
    
    % get the distance from left end of the shaft to each node 
    nodeDistance{iShaft} = [];
    for iSegment = 1:1:rowSegmentNum(iShaft)
        nodesInSegment = linspace( keyPoints{iShaft}(iSegment),...
                                   keyPoints{iShaft}(iSegment+1),...
                                   elementNum(iSegment)+1 );
        nodeDistance{iShaft} = [nodeDistance{iShaft}, nodesInSegment];
    end
    
    
    % merge the same value
    nodeDistance{iShaft} = unique(nodeDistance{iShaft});
    
end % end for iShaft = 1:1:Shaft.amount

%%

% generate a new struct saving the information of each node: Node
nodeNum = sum( cellfun(@length,nodeDistance) );
Node = struct(  'name',             cell(nodeNum,1),...
                'onShaftNo',        [],...
                'onShaftDistance',  [],...
                'diskNo',           [],...
                'bearingNo',        [],...
                'isBearing',        [],...
                'dof',              [] ); %initial

for iNode = 1:1:nodeNum
    if hasRub
        [Node(iNode).rubImpactNo] = [];
    end

    if hasInterBearing
        Node(iNode).interBearingNo = [];
    end

    if hasCoupling
        Node(iNode).couplingNo = [];
    end
    
    if hasLoosingBearing
        Node(iNode).isLoosingBearing = [];
    end
end

 
% check every node            
iShaft = 1;
previousShaftNodeNum = 0;
for iNode = 1:1:nodeNum
    Node(iNode).name = iNode;
    
    % Node.onShaftNo and Node.onShaftDistance
    Node(iNode).onShaftNo = iShaft;
    distanceHere = nodeDistance{iShaft}(iNode - previousShaftNodeNum);
    Node(iNode).onShaftDistance = distanceHere;
    Node(iNode).dof = Shaft.dofOfEachNodes(iShaft);
    
    % judge if the shaft is end
    [previousShaftNodeNum, iShaft] = judgeShaftEnd(previousShaftNodeNum,...
                                     nodeDistance, iShaft, iNode); 
end % end for iNode

%%

% Node.diskNo and Disk.positionOnShaftNode
[Disk,Node] = matchElement(Disk, 'disk', Node, nodeDistance,nodeNum,...
              THRESHOLD_COEFFICIENT, Shaft);
           
%%

% Node.bearingNo and Bearing.positionOnShaftNode
[Bearing,Node] = matchElement(Bearing, 'bearing', Node, nodeDistance,nodeNum,...
                 THRESHOLD_COEFFICIENT, Shaft);

%%

% Node.rubImpactNo and RubImpact.positionOnShaftNode
if hasRub
    [RubImpact,Node] = matchElement(RubImpact, 'rubImpact', Node,...
                       nodeDistance,nodeNum,THRESHOLD_COEFFICIENT, Shaft);
end

%%

% Node.couplingNo and Coupling.positionOnShaftNode
if hasCoupling
    [Coupling,Node] = matchElement(Coupling, 'couplingMisalignment', Node,...
                      nodeDistance,nodeNum,THRESHOLD_COEFFICIENT, Shaft);
end

%%

% Node.interBearingNo and InterBearing.positionOnShaftNode
if hasInterBearing
    [InterBearing,Node] = matchElement(InterBearing, 'interBearing', Node,...
                          nodeDistance,nodeNum,THRESHOLD_COEFFICIENT, Shaft);
end

%%

% Node.isLoosingBearing
if hasLoosingBearing
    loosingBearingIndex = LoosingBearing.inBearingNo;
    loosingNode = Bearing.positionOnShaftNode(loosingBearingIndex);
    Node(loosingNode).isLoosingBearing = deal(true);
    for iNode = 1:1:nodeNum
        if isempty(Node(iNode).isLoosingBearing)
            Node(iNode).isLoosingBearing = false;
        end
    end
end


%%

% Node.isBearing and Bearing.positionNode
% generater bearing nodes
isInter = false;
[nodeNum, Bearing, Node] = addNode(nodeNum, Bearing, Node, hasLoosingBearing, false);

%%

% Node.isInterBearing and InterBearing.positionNode
% generater bearing nodes
isInter = true;
[nodeNum, InterBearing, Node] = addNode(nodeNum, InterBearing, Node, hasLoosingBearing, isInter);

%%

% Mesh.dofInterval
dofOfEachNode = [Node.dof];
dofInterval = zeros(nodeNum,2);
for iNode = 1:1:nodeNum
    endDof = sum( dofOfEachNode(1:iNode) );
    startDof = endDof - dofOfEachNode(iNode) +1;
    dofInterval(iNode,:) = [startDof, endDof];
end

%%

% Mesh.dofOnNodeNo
dofOnNodeNo = zeros(sum([Node.dof]),1);
dofNo = 1;
for iNode = 1:1:nodeNum
    for iDof = 1:1:Node(iNode).dof
        dofOnNodeNo(dofNo) = Node(iNode).name;
        dofNo = dofNo + 1;
    end
end

%% 

% output
Parameter = InitialParameter;


% output struct: Mesh
Mesh.Node = Node;
Mesh.keyPointsDistance = keyPoints;
Mesh.nodeDistance = nodeDistance;
Mesh.nodeNum = nodeNum;
Mesh.dofNum = sum([Node.dof]);
Mesh.dofOnNodeNo = dofOnNodeNo;
Mesh.dofInterval = dofInterval;
Parameter.Mesh = Mesh;


%output the node position of every important elements
Parameter.Disk = Disk;
Parameter.Bearing = Bearing;

if hasRub
    Parameter.RubImpact = RubImpact;
end

if hasCoupling
    Parameter.CouplingMisalignment = Coupling;
end

if hasInterBearing
    Parameter.IntermediateBearing = InterBearing;
end


%%

% subfunction 1
% write the node information into important elements (e.g. Disk, Bearing)
% and the struct Node
function [Element,Node] = matchElement(Element, elementName, Node,...
                          nodeDistance, nodeNum, THRESHOLD_COEFFICIENT, Shaft)
            
if strcmp(elementName,'interBearing')
    columnNum = 2;
else
    columnNum = 1;
end

Element.positionOnShaftNode = zeros(Element.amount,columnNum);    
shaftNo = 1;
previousShaftNode = 0;
iElement = 1;
for nodeNo = 1:1:nodeNum
    % calculate the distance at nodeNo
    distanceH = nodeDistance{shaftNo}(nodeNo - previousShaftNode);
    
    
    % judge if this node has important element (e.g. disk, bearing)
    if strcmp(elementName,'interBearing')
        isShaftHere = shaftNo == Element.betweenShaftNo;
    else
        isShaftHere = shaftNo == Element.inShaftNo;
    end
    isDistanceHere = abs(distanceH - Element.positionOnShaftDistance)...
                     < Shaft.totalLength(shaftNo)/THRESHOLD_COEFFICIENT;
    isElementHere = isShaftHere & isDistanceHere; % logical matrix
    isElementHereNum = length(isElementHere(isElementHere == true));

    
    % write the information into Element and Node
    if isElementHereNum == 1 % there is a single disk on this position
        
        % Element.positionOnShaftNode 
        indexElement = find(isElementHere == true);
        Element.positionOnShaftNode(indexElement) = nodeNo;
        
        % Node.xxxxNo
        switch elementName
            case 'disk'
                Node(nodeNo).diskNo = iElement;
            case 'bearing'
                Node(nodeNo).bearingNo = iElement;
            case 'rubImpact'
                Node(nodeNo).rubImpactNo = iElement;
            case 'couplingMisalignment'
                Node(nodeNo).couplingNo = iElement;
            case 'interBearing'
                % only write the interBearing No. of 1st column of interBearing
                if indexElement(:,1) == 1 
                    Node(nodeNo).interBearingNo = iElement; 
                else
                    iElement = iElement - 1;
                end % end if
            otherwise
                error('wrong elementName in matchElement()')
        end % end switch
        
        iElement = iElement +1;
        
    elseif isElementHereNum > 1 % there are too many disks on this position
       
        error(['the distance between disks is too close, please modify the',...
               'THRESHOLD_COEFFICIENT in meshModel() or redefine the position of disks']);

    end % end if isElementHereNum == 1
    
    
    % judge if the shaft is end
    [previousShaftNode, shaftNo] = judgeShaftEnd(previousShaftNode,...
                                   nodeDistance, shaftNo, nodeNo);

end % end for 


% write the interBearing No. of 2nd column of interBearing
if strcmp(elementName,'interBearing')
        index1Column = Element.positionOnShaftNode(:,1);
        index2Column = Element.positionOnShaftNode(:,2);
        Node(index2Column).interBearingNo = Node(index1Column).interBearingNo;
end % end if

end % end subfunction matchElement()

%%

% subfunction 2
% judge if this node is the end of the shaft. true -> iShift +1 and
% previousShaftNodeNum + the number of nodes last shaft having
function [previousShaftNodeNum, iShaft]...
        = judgeShaftEnd(previousShaftNodeNum, nodeDistance, iShaft, iNode)

isShaftEnd = iNode == (previousShaftNodeNum +...
                   length(nodeDistance{iShaft}) );
if isShaftEnd
    previousShaftNodeNum = previousShaftNodeNum +...
                           length(nodeDistance{iShaft});
    iShaft = iShaft +1;
end % end if

end % end sub function judgeShaftEnd()

%%

% subfunction 3
% add Node.isInterBearing and InterBearing.positionNode, Node.isBearing and
% Bearing.positionNode for Bearing and InterBearing
function [nodeNum, Element, Node] = addNode(nodeNum, Element, Node, hasLoosingBearing, isInter)
indexHasMass = find(Element.mass' ~= 0); % index in 1-d
[indexHasMassRow, indexHasMassCol] = find(Element.mass' ~= 0); % index in 2-d
bearingNodeNum = length(indexHasMass);
bearingMassColumnNum = size(Element.mass,2);
Element.positionNode = zeros(Element.amount,bearingMassColumnNum); % initial
for iBearingNode = 1:1:bearingNodeNum
    newNodeNo = nodeNum + iBearingNode;
    onNodeNo = Element.positionOnShaftNode(indexHasMassCol(iBearingNode),:);
    % copy information in Node
    Node(newNodeNo).onShaftNo = [Node(onNodeNo).onShaftNo];
    Node(newNodeNo).onShaftDistance = [Node(onNodeNo).onShaftDistance];
    
    if isInter
        Node(newNodeNo).interBearingNo = unique([Node(onNodeNo).interBearingNo]);
    else
        Node(newNodeNo).bearingNo = Node(onNodeNo).bearingNo;
    end
    
    if hasLoosingBearing
        Node(newNodeNo).isLoosingBearing = Node(onNodeNo).isLoosingBearing;
    end
    % write new information in Node
    Node(newNodeNo).name = newNodeNo; 
    Node(newNodeNo).dof = Element.dofOfEachNodes(indexHasMassCol(iBearingNode),indexHasMassRow(iBearingNode));
    Node(newNodeNo).isBearing = true;
    % Element.positionNode
    Element.positionNode(indexHasMassCol(iBearingNode),indexHasMassRow(iBearingNode)) = newNodeNo;
end


for iNodee = 1:1:nodeNum
    if isempty(Node(iNodee).isBearing)
        Node(iNodee).isBearing = false;
    end
end
nodeNum = nodeNum + bearingNodeNum;
    
end
end % end main function