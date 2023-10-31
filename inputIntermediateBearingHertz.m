%% inputIntermediateBearing
% Input parameters about intermediate bearing
%% Syntax
%  OutputParameter = InputIntermediateBearing(InputParameter)
%% Description
% InputParameter is a struct data saving all initial parameters except
% intermediate bearing parameters.
% 
% OutputParameter is a strct data recording all InputParameter and the data
% about intermediate bearing.
%
% All of parameters about intermediate bearing should be typed into this .m file.


%%
function OutputParameter = inputIntermediateBearingHertz(InputParameter)

% typing the parameter about intermediate bearing
IntermediateBearing.amount          = 1;
% shaft no. connected by same bearing in row; different bearings in column
IntermediateBearing.betweenShaftNo  = [2 1]; % n*2
% the same bearing in row; different bearings in column, n*2
IntermediateBearing.positionOnShaftDistance = 1e-3 * [0 258];
IntermediateBearing.stiffness       = [1e9]; % in column, n*1
IntermediateBearing.damping         = [300]; % in column, n*1

%%

%check the input data
if size(IntermediateBearing.betweenShaftNo, 2) ~= 2
    error('IntermediateBearing.betweenShaftNo must be a 2 column matrix')
end
               
if size(IntermediateBearing.betweenShaftNo, 2) ~= 2
    error('IntermediateBearing.positionOnShaftDistance must be a 2 column matrix')
end

checkInputData(IntermediateBearing);

%%

% sort the betweenShaftNo and PositionOnShaftDistance (each row)
shaftNo = IntermediateBearing.betweenShaftNo; % short the variable
position = IntermediateBearing.positionOnShaftDistance;
for iBearing = 1:1:IntermediateBearing.amount
    if shaftNo(iBearing,1)>shaftNo(iBearing,2)
       % exchange the column 1 and column 2 in iBearing row 
       temporary = shaftNo(iBearing,1);
       shaftNo(iBearing,1) = shaftNo(iBearing,2);
       shaftNo(iBearing,2) = temporary;
       temporary = position(iBearing,1);
       position(iBearing,1) = position(iBearing,2);
       position(iBearing,2) = temporary;
    end
end

% sort the betweenShaftNo and PositionOnShaftDistance (each column)
[shaftNo,index] = sortrows(shaftNo,1);
position = position(index,:);
stiffness = IntermediateBearing.stiffness(index,:);
damping = IntermediateBearing.damping(index,:);

IntermediateBearing.betweenShaftNo = shaftNo;
IntermediateBearing.positionOnShaftDistance = position;
IntermediateBearing.stiffness = stiffness;
IntermediateBearing.damping = damping;
%%

OutputParameter = InputParameter;
OutputParameter.IntermediateBearing = IntermediateBearing;
OutputParameter.ComponentSwitch.hasIntermediateBearing = true;
end % for function inputIntermediateBearing()




