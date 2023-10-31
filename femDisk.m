%% femDisk
% generate the globe mass, stiffness, gyroscopic matrix of disks
%% Syntax
% [M, G, N] = femDisk(Disk,nodeDof)
%% Description
% Disk is a struct saving the physical parameters of Disks with fields:
% amount, dofOfEachNodes, radius, thickness, density, positionOnShaftNode
%
% nodeDof: is a array (the number of nodes  * 1) saving the dof of each 
% node.
%
% M, K, G, N are mass, stiffness, gyroscopic, N matrix of disks. (n*n,
% n is the number of all nodes)



function [M, G, N, Q] = femDisk(Disk,nodeDof)


% generate elements
Me = cell(Disk.amount,1); 
Ge = cell(Disk.amount,1);
Ne = cell(Disk.amount,1);
Temporary = rmfield(Disk,'amount'); % for extract part of data of Shaft

for iDisk = 1:1:Disk.amount
    % get the information of ith Disk
    ADisk = getStructPiece(Temporary,iDisk,[]);
    % generate elements
    [Me{iDisk}, Ge{iDisk}, Ne{iDisk}] = diskElement(ADisk); 
end

%%

% generate global matrices
dofNum = sum(nodeDof);
M = zeros(dofNum, dofNum);
G = zeros(dofNum, dofNum);
N = zeros(dofNum, dofNum);

%%

% calculate the position of disk element in global matrix
diskOnDofPosition = findIndex(Disk.positionOnShaftNode,nodeDof);

%%

% put disk elements into global matrices
for iDisk = 1:1:Disk.amount
   M = addElementIn( M, Me{iDisk}, diskOnDofPosition(iDisk, :) );
   G = addElementIn( G, Ge{iDisk}, diskOnDofPosition(iDisk, :) );
   N = addElementIn( N, Ne{iDisk}, diskOnDofPosition(iDisk, :) );
end

Q = zeros(dofNum,1);

end