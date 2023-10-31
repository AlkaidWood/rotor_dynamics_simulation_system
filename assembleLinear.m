%% assembleLinear
% generate the globe mass, stiffness, gyroscopic matrix of shafts
%% Syntax
% B = assembleLinear(A, intersectionRow, intersectionColumn)
%% Description
% A is a cell (n*1). Each cell saves different 2_D matrices.
%
% intersectionRow/Column is a integer vector, indicates the dimension the 
% matrices intersecting.
%
% B is a matrix after intersecting

function B = assembleLinear(A, intersectionRow, intersectionColumn)

% check input data
matrixNum = length(A);
indexNumRow = length(intersectionRow);
indexNumColumn = length(intersectionColumn);
if ~isequal(matrixNum, indexNumRow+1, indexNumColumn+1)
    error('the dimensions of input parameters are different');
end % end if 

%%

% calculate the size of matrix B and initial
[rowNumA, columnNumA] = cellfun(@size, A);
rowNumB = sum(rowNumA) - sum(intersectionRow);
columnNumB = sum(columnNumA) - sum(intersectionColumn);
B = zeros(rowNumB, columnNumB);

%%

% intersecte matrices
for iMatrix = 1:1:matrixNum
    [rowMatrix, columnMatrix] = size(A{iMatrix});
    rowIndexA = 1:1:rowMatrix;
    columnIndexA = 1:1:columnMatrix;
    if iMatrix == 1
        B(rowIndexA,columnIndexA) = A{iMatrix}(rowIndexA,columnIndexA);
        [endRow, endColumn] = size(A{iMatrix});
    else
        positionRow = endRow - intersectionRow(iMatrix-1);
        positionColumn = endColumn - intersectionColumn(iMatrix-1);
        rowIndexB = rowIndexA + positionRow;
        columnIndexB = columnIndexA + positionColumn;
        B(rowIndexB,columnIndexB) = B(rowIndexB,columnIndexB)...
                                    + A{iMatrix}(rowIndexA,columnIndexA);
        endRow = positionRow + rowMatrix;
        endColumn = positionColumn + columnMatrix;
    end % end if
end % end for

end % end function