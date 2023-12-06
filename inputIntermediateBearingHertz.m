%% inputIntermediateBearingHertz
% Input parameters about intermediate bearing considering Hertzian contact
%% Syntax
% OutputParameter = inputIntermediateBearingHertz(InputParameter)
%% Description
% InputParameter is a struct data saving all initial parameters except
% intermediate bearing parameters considering Hertzian contact.
% 
% OutputParameter is a strct data recording all InputParameter and the data
% about intermediate bearing considering Hertzian contact.
%
% All of parameters about intermediate bearing considering Hertzian contact
% should be typed into this .m file.
%
% All vectors in this .m file should be in column


%%
function OutputParameter = inputIntermediateBearingHertz(InputParameter)

% typing the parameter about intermediate bearing
IntermediateBearingH.amount          = 1;
% shaft no. connected by same bearing in row; different bearings in column
IntermediateBearingH.betweenShaftNo  = [2 1]; % n*2
% dof
IntermediateBearingH.dofOfEachNodes = [2];% if mass=0, dof must be 0 
% the same bearing in row; different bearings in column, n*2
IntermediateBearingH.positionOnShaftDistance = 1e-3 * [0 258];
IntermediateBearingH.isHertzian      = [true]; % boolean
% M K C, elements in the same row: the MKC at the same position of the
% shaft; mass(1,1) -> mass(1,n): the mass near the outer shaft -> the mass 
% near the inner shaft
IntermediateBearingH.stiffness       = [1e9]; % N/m, in column, n*1
IntermediateBearingH.damping         = [300]; % N/s^2, in column, n*1
IntermediateBearingH.mass            = [3]; % kg
IntermediateBearingH.rollerNum        = [8];
IntermediateBearingH.radiusInnerRace = [28.8e-3]; % m
IntermediateBearingH.radiusOuterRace = [47e-3]; % m
IntermediateBearingH.clearance = [10e-6]; % m
IntermediateBearingH.contactStiffness = [1.08e10]; % N*m^-3/2
IntermediateBearingH.coefficient = [3/2]; % =3/2 in a ball bearing; = 10/9 in a roller bearing

%%

%check the input data
if size(IntermediateBearingH.betweenShaftNo, 2) ~= 2
    error('IntermediateBearingH.betweenShaftNo must be a 2 column matrix')
end
               
if size(IntermediateBearingH.positionOnShaftDistance, 2) ~= 2
    error('IntermediateBearingH.positionOnShaftDistance must be a 2 column matrix')
end

checkInputData(IntermediateBearingH);

%%

% sort the betweenShaftNo and PositionOnShaftDistance (sort column)
shaftNo = IntermediateBearingH.betweenShaftNo; % short the variable
position = IntermediateBearingH.positionOnShaftDistance;
for iBearing = 1:1:IntermediateBearingH.amount
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

% sort other parameters with the shaftNo (sort row)
[shaftNo,index] = sortrows(shaftNo,1);
position = position(index,:);
stiffness = IntermediateBearingH.stiffness(index,:);
damping = IntermediateBearingH.damping(index,:);

IntermediateBearingH.betweenShaftNo = shaftNo;
IntermediateBearingH.positionOnShaftDistance = position;
IntermediateBearingH.stiffness = stiffness;
IntermediateBearingH.damping = damping;
%%

OutputParameter = InputParameter;
OutputParameter.IntermediateBearing = IntermediateBearingH;
OutputParameter.ComponentSwitch.hasIntermediateBearing = true;
end % for function inputIntermediateBearingHertz()




