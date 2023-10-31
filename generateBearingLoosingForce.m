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
for iNode = 1:1:Mesh.nodeNum
    isBearing = Node(iNode).isBearing == true;
    isLoosingBearing = Node(iNode).bearingNo == LoosingBearing.inBearingNo;
    if isBearing && isLoosingBearing
        loosingNodeNo = Node(iNode).name;
    end
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