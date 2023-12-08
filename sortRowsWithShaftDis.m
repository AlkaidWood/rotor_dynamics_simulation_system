%% sortRowsWithShaftDis
% Input Bearing or IntermediateBearing (struct in initial parameters)
%% Syntax
%  orderedElement = sortRowsWithShaftDis(Element)
%% Description
% Element: Bearing or IntermediateBearing (struct in initial parameters)
% 
% orderedElement: return a struct with order of both shaft and displacement
% 
% e.g. Bearing.onShaftNo = [1; 2; 1]; Bearing.positionOnShaftDistance = [5;
% 2; 1]; Output => Bearing.onShaftNo = [1; 1; 2];
% Bearing.positionOnShaftDistance = [1; 5; 2];

function orderedElement = sortRowsWithShaftDis(Element)
Temporary = rmfield(Element,'amount');

% divide Element into pieces respect to shaft No.
if isfield(Temporary, 'inShaftNo')
    inShaftNo = Temporary.inShaftNo(:,1);
else
    inShaftNo = Temporary.betweenShaftNo(:,1);
end
shaftNoSet = unique(inShaftNo);
shaftNum = length(shaftNoSet);
% convert strct to cell
fieldName = fieldnames(Temporary); % save the field names
cellData = struct2cell(Temporary); % save the data
% find the index of certain field name in the cell data
if isfield(Temporary, 'inShaftNo')
    inShaftNoIndex = ismember(fieldName, 'inShaftNo');
else
    inShaftNoIndex = ismember(fieldName, 'betweenShaftNo');
end
distanceIndex = ismember(fieldName, 'positionOnShaftDistance');
orderCellData = cell(length(cellData),1);
for iShaft = 1:1:shaftNum
    % find columns with same shaft no
    index = cellData{inShaftNoIndex}(:,1)==iShaft;
    % extract the data
    cellDataShaft = cell(length(cellData),1);
    for iRow = 1:1:length(cellData)
        cellDataShaft{iRow} = cellData{iRow}(index,:);
    end
    % order the data respect to distance
    [~,indexDis] = sort(cellDataShaft{distanceIndex}(:,1));
    for iRow = 1:1:length(cellData)
        cellDataShaft{iRow} = cellDataShaft{iRow}(indexDis,:);
        % save the ordered data
        orderCellData{iRow} = [orderCellData{iRow};cellDataShaft{iRow}];
    end % end for iRow
end % end for iShaft

% convert cell to struct
Temporary = cell2struct(orderCellData, fieldName, 1);

% output
orderedElement = Temporary;
orderedElement.amount = Element.amount;

end