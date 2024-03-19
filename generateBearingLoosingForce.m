%% generateBearingLoosingForce
% generate a 'bearingLoosingForce.m' file
%% Syntax
% generateBearingLoosingForce(LoosingBearing)
%% Description
% LoosingBearing: is a struct saving the relevant data
%
% generate a bearingLoosingForce.m file in the root directory


function generateBearingLoosingForce(LoosingBearing, Mesh)

% check the exist of bearingLoosingForce.m and create .txt
if isfile('bearingLoosingForce.m')
    delete bearingLoosingForce.m
end

blf = fopen('bearingLoosingForce.txt','w');

%%

% write comments, firstly

comments = [...
"%% bearingLoosingForce";...
"% saving the equation of loosing bearing force in this function";...
"%% Syntax";...
"% [K, C] = bearingLoosingForce(qn)";...
"%% Description";...
"% qn: is the displacement at the n-th time";...
"% ";...
"% Matrix: is a struct saving all matrix used in dynamic equation";...
"% ";...
"% K, C: are matrix after loosing";...
" ";...
" "...
];

%%

% function start
functionStart = [...
"function [K, C] = bearingLoosingForce(qn, Matrix)";...
" "...
];

fprintf(blf,'%s\n',comments);
fprintf(blf,'%s\n',functionStart);

%%

% calculate the loosing dof
Node = Mesh.Node;
loosingNodeNo = zeros(1000,1); % initial
loosingNum = 0; % record the mass number of the loosing bearing

% If the LoosingBearing.loosingPositionNo == 1:
% k1,c1 will be change to loosing k,c when the y-displacement of the shaft
% connecting the loosing bearing exceeds the interval
% If the LoosingBearing.loosingPositionNo = j, j > 1:
% kj,cj will be change to loosing k,c when the y-displacement of the
% (j-1)th mass of the loosing bearing exceeds the interval
for iNode = 1:1:Mesh.nodeNum
    if LoosingBearing.loosingPositionNo == 1 % denote that the k1,c1 is loosing
        isElement = Node(iNode).isBearing == false; % this element is shaft
    else % denote the kj,cj is loosing, j>1
        isElement = Node(iNode).isBearing == true; % this element is bearing
    end % end if
    isLoosingBearing = Node(iNode).isLoosingBearing;
    if isElement && isLoosingBearing
        loosingNum = loosingNum + 1;
        loosingNodeNo(loosingNum) = Node(iNode).name;
    end % end if
end % end for
loosingNodeNo = loosingNodeNo(loosingNodeNo~=0); % delete the 0 elements
if LoosingBearing.loosingPositionNo ~= 1
    loosingNodeNo = loosingNodeNo(LoosingBearing.loosingPositionNo-1);
end
loosingDof = Mesh.dofInterval(loosingNodeNo,1) + 1; % loosing direction on y
interval = LoosingBearing.interval;

%%

% function body
functionBody = {...
 ' ';...
['if qn(', num2str(loosingDof), ') >= 0 && qn(', num2str(loosingDof), ') <= ', num2str(interval)];...
 '    K = Matrix.stiffnessLoosing;';...
 '    C = Matrix.dampingLoosing;';...
 'else';...
 '    K = Matrix.stiffness;';...
 '    C = Matrix.damping;';...
 'end';...
 ' ';...
};

functionBody = cell2string(functionBody);

fprintf(blf,'%s\n', functionBody);

%%

% function end
functionEnd = [...
"end";...
" "...
];
fprintf(blf,'%s\n',functionEnd);

%%

% close .txt and transfer .txt -> .m
fclose(blf);
system('rename bearingLoosingForce.txt bearingLoosingForce.m');
end