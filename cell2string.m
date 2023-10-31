
function stringData = cell2string(cellData)
stringData = zeros(size(cellData));
stringData = string(stringData);

for iRow = 1:1:size(cellData,1)
    for iColumn =1:1:size(cellData,2)
        stringData(iRow,iColumn) = cellData{iRow,iColumn};
    end % end for iColumn
end % end for iRow

end  % end sub function