%% femInterBearing
% generate the globe mass, stiffness, damping matrix of bearings
%% Syntax
% [K, C] = femInterBearing( InterBearing, nodeDof )
%% Description
% InterBearing is a struct saving the physical parameters of bearings with 
% fields: amount, stiffness, damping, positionOnShaftNode
%
% nodeDof: is a array (the number of nodes  * 1) saving the dof of each 
% node.
%
% K, C are stiffness, damping matrix of intermediate bearings. (n*n,
% n is the number of all nodes)

function  [K, C] = femInterBearing( InterBearing, nodeDof )

% generate global matrices
dofNum = sum(nodeDof);
K = zeros(dofNum, dofNum);
C = zeros(dofNum, dofNum);

%%
% gnerate elements
Ke = cell(InterBearing.amount,1); 
Ce = cell(InterBearing.amount,1);
Temporary = rmfield(InterBearing,'amount');

for iBearing = 1:1:InterBearing.amount
    % get the information of ith intermediate braring
    ABearing = getStructPiece(Temporary,iBearing,[],false); % a normal bearing
    ABearing.dofOnShaftNode = nodeDof(InterBearing.positionOnShaftNode(iBearing,:));
    % generate elements 
    [Ke{iBearing}, Ce{iBearing}] = bearingElementInter(ABearing); 
end

%%

% find index of elements in the global matrix
bearingIndex = findIndex(InterBearing.positionOnShaftNode,nodeDof); 

%%

% put the intermediate bearing elements into global matrix
for iBearing = 1:1:InterBearing.amount
    K = repeatAdd(K, Ke, iBearing, bearingIndex);
    C = repeatAdd(C, Ce, iBearing, bearingIndex);
end

%%

% sub function
function B = repeatAdd(A, Ae, iObject, aIndex)
    A = addElementIn(A, Ae{iObject}{1,1}, aIndex(2*iObject-1,[1,2]));
    A = addElementIn(A, Ae{iObject}{1,2}, aIndex(2*iObject-1,[3,4]));
    A = addElementIn(A, Ae{iObject}{2,1}, aIndex(2*iObject,[1,2]));
    A = addElementIn(A, Ae{iObject}{2,2}, aIndex(2*iObject,[3,4]));
    B = A;
end % end sub function

end