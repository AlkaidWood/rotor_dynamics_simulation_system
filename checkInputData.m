%% checkInputData
% Estimating the dimension of the inputdata, providing warning
%% Syntax
% checkInputData(checkData)
%% Description
% checkData is a struct data saving a part of initial parameters. 
% checkData must have a fieldname: amount 
%
% The column dimension of other fields in checkData should equal
% checkData.amount, otherwise the function outputs an error information.


function checkInputData(checkData)

rowNum = checkData.amount;
checkData = rmfield(checkData,'amount'); % amount is scalar, no check
checkData = struct2cell(checkData);

for iData = 1:1:size(checkData,1)
    if size(checkData{iData}, 1) ~= rowNum 
        % get first input variable name of this function
        variableName = inputname(1); 
        stackName = dbstack(1); % get stack information
        % get the function name calling this function
        fileName = stackName(1).name; 
        error([ 'In function ', fileName  '(), ', 'the input data in <', variableName, '> must be a ', num2str(rowNum), ' row matrix.'])
    end % for if
end % for loop

end % for subfunction checkInputData()