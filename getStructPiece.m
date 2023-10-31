%% getStructPiece
% get a piece of data of target struct (the dimension of data in each field
% of target struct must be same)
%% Syntax
% StructPiece = getStructPiece(Target,rowIndex,columnIndex)
%% Description
% Target: is a struct whose data will be extracted
%
% rowIndex/column: is a int indicating the row/column index of data in
% every field which will be extrcted
%
% StructPiece: is a struct including the same fields with Target and having
% a piece of data of Target.


function StructPiece = getStructPiece(Target,rowIndex,columnIndex, isCheck)

if nargin < 4
    isCheck = true;
end

%%

% check the input
checkData = struct2cell(Target);
fieldNum = length(checkData);
dimension = zeros(fieldNum,1);

for iData =1:1:fieldNum
    dimension(iData) = length( checkData{iData} );
end

isDimensionEqual = length(unique(dimension)) == 1;
if ~isDimensionEqual && isCheck
    error('the dimension of data in every field must be same')
end

isNullRow = isempty(rowIndex);
isNullColumn = isempty(columnIndex);

if isNullRow  && isNullColumn
    error('index is empty')
end



%%

fieldName = fieldnames(Target);

for iField = 1:1:fieldNum
    fullData = Target.(fieldName{iField});
    if isNullRow
        StructPiece.(fieldName{iField}) = fullData(:,columnIndex);
    elseif isNullColumn
        StructPiece.(fieldName{iField}) = fullData(rowIndex,:);
    else
        StructPiece.(fieldName{iField}) = fullData(rowIndex,columnIndex);
    end % end if 
end % end for

end % end function