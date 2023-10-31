%% addElementIn
% add a matrix into another lager matrix
%% Syntax
% C = addElementIn(A, B, position)
%% Description
% A is the a matrix lagger (dimension) than B .
%
% B is the matrix waiting to be added
%
% positiion = [x,y] is a index where the first elements of B should be added
%
% C is the matrix after combining A and B

function C = addElementIn(A, B, position)

% check input
[rowNumA, columnNumA] = size(A);
[rowNumB, columnNumB] = size(B);
isOutRow = rowNumB+position(1)-1 > rowNumA;
isOutColumn = columnNumB+position(2)-1 > columnNumA;
if isOutRow || isOutColumn
   error('the smaller matrix exceed the boundary of lagger matrix'); 
end
if (length(position) ~= 2) || (length(unique( size(position) )) ~= 2)
    error('the position should be a 1*2 or 2*1 array')
end

%%
% add element of B into A
rowIndexInB = 1:1:rowNumB;
columnIndexInB = 1:1:columnNumB;
rowIndexInA = rowIndexInB + position(1) - 1;
columnIndexInA = columnIndexInB + position(2) -1;


C = A;
C(rowIndexInA, columnIndexInA) = A(rowIndexInA, columnIndexInA)...
                                 + B(rowIndexInB, columnIndexInB);
                             
end