%% findIndex
% find the index of element in global matrix
%% Syntax
% elementIndex = findIndex(positionNode,nodeDof)
%% Description
% positionNode: is a column or n*2 matrix saving the nodes where elements 
% are. 
% 
% (1)positionNode is column: the nodes where elements are. 
% 
% (2)positionNode is n*2 matrix: this element 
% contains two nonadjacent nodes; the first column contains the first nodes 
% and the second column contains the second nodes (for bearing-mass 
% elemeents, intermediate elements; there are four matrix waiting for
% putting into global matrix.
%
% nodeDof: is a array (the number of nodes  * 1) saving the dof of each 
% node.
%
% elementIndex: (1) is a n*2 matrix where each row indicates the index of
% the element in the global matrix (when positionNode is a column);
% (2) is a 2n * 4 matrix contains 4*n index:
%
% elementIndex = [ index11, index12;...    element 1
%                  index21, index22;...    element 1
%                  index11, index12;...    element 2
%                  index21, index22 ];     element 2

function elementIndex = findIndex(positionNode,nodeDof)

% check input
if size(positionNode,2)>2
    error('incorrect dimension of first input parameter');
end

%%
% calculate the position of element in global matrix
elementNum = size(positionNode,1);

if size(positionNode,2) == 1
    
    node1No = positionNode;
    elementIndex = zeros(elementNum, 2);
    for iElement = 1:1:elementNum
        num1 = sum(nodeDof(1:node1No(iElement)-1)) + 1;
        % position of the first disk element located; 
        elementIndex(iElement,1) = num1;
        elementIndex(iElement,2) = num1;
    end % end for
    
elseif size(positionNode,2) == 2
    
    node1No = positionNode(:,1);
    node2No = positionNode(:,2);
    elementIndex = zeros(2*elementNum,4);
    for iElement = 1:1:elementNum
        num1 = sum(nodeDof(1:node1No(iElement)-1)) + 1;
        num2 = sum(nodeDof(1:node2No(iElement)-1)) + 1;
        % Index11
        elementIndex(2*iElement-1,1) = num1;
        elementIndex(2*iElement-1,2) = num1;
        % Index12
        elementIndex(2*iElement-1,3) = num1;
        elementIndex(2*iElement-1,4) = num2;
        % Index21
        elementIndex(2*iElement,1) = num2;
        elementIndex(2*iElement,2) = num1;
        % Index22
        elementIndex(2*iElement,3) = num2;
        elementIndex(2*iElement,4) = num2;
    end % end for  
    
else
    
    error('incorrect dimension of first input parameter');
    
end % end if


end % end funcion

