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
IntermediateBearingH.amount          = 4;
% shaft no. connected by same bearing in row; different bearings in column
IntermediateBearingH.betweenShaftNo  = [2, 1;...
                                        1, 2;...
                                        1, 2;...
                                        1, 2]; % n*2
% dof
IntermediateBearingH.dofOfEachNodes = [2,0,0;...
                                       2,2,0;...
                                       2,2,2;...
                                       0,0,0];% if mass=0, dof must be 0 
% the same bearing in row; different bearings in column, n*2
IntermediateBearingH.positionOnShaftDistance = 1e-3 * [0,   258;...
                                                       400, 50;...
                                                       500, 100;...
                                                       600, 300];
IntermediateBearingH.isHertzian      = [true;false;false;true]; % boolean
% M K C, elements in the same row: the MKC at the same position of the
% shaft; mass(1,1) -> mass(1,n): 
% the mass near the inner shaft -》the mass near the outer shaft
% If isHertizian and no mass, the corresponding k c will be added in
% global matrix normally; the model:
% Inner shaft--Hertz+k1c1--Outer shaft;
% If is no Hertzian and with mass: there are n mass in a row, and n+1 k c 
% for a bearing; the model will be established as:
% Inner shaft--k1c1--m1--k2c2--m2--k3c3--m3--k4c4- ...-mn--k(n+1)c(n+1)--Outer shaft;
% If is no Hertzian and no mass, the k c will be added in global 
% matrix normally; the model:
% Inner shaft--k1c1--Outer shaft;
% If isHertzian and with mass, the hertzian force will be added at the mass
% in the first column (near the shaft); the model:
% Inner shaft--Hertz+k1c1--m1--k2c2--m2--k3c3--m3--k4c4- ...-mn--k(n+1)c(n+1)--Outer shaft;
IntermediateBearingH.stiffness       = [0,   1e9, 0,   0;...
                                        1e3, 1e3, 1e3, 0;...
                                        1e3, 1e3, 1e3, 1e3;...
                                        1e3, 0,   0,   0]; % N/m, in column, n*1
IntermediateBearingH.damping         = [0,   300, 0,   0;...
                                        50,  100, 150, 0;...
                                        50,  100, 150, 300;...
                                        300, 0,   0,   0]; % N/s^2, in column, n*1
IntermediateBearingH.mass            = [3, 0, 0;...
                                        3, 4, 0;...
                                        3, 4, 5;...
                                        0, 0, 0]; % kg
% if there is no Hertizan contact force, please set n*1 zero vector for following parameters, where n is the number of the intermediate bearing                                 
IntermediateBearingH.rollerNum        = [8;0;0;8];
IntermediateBearingH.radiusInnerRace = [28.8e-3;0;0;28.8e-3]; % m
IntermediateBearingH.radiusOuterRace = [47e-3;0;0;47e-3]; % m
IntermediateBearingH.clearance = [10e-6;0;0;10e-6]; % m
IntermediateBearingH.contactStiffness = [1.08e10;0;0;1.08e10]; % N*m^-3/2
IntermediateBearingH.coefficient = [3/2;0;0;3/2]; % =3/2 in a ball bearing; = 10/9 in a roller bearing

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
IntermediateBearingH.betweenShaftNo = shaftNo;
IntermediateBearingH.positionOnShaftDistance = position;
% sort columns in struct
IntermediateBearingH = sortRowsWithShaftDis(IntermediateBearingH);

%%

OutputParameter = InputParameter;
OutputParameter.IntermediateBearing = IntermediateBearingH;
OutputParameter.ComponentSwitch.hasIntermediateBearing = true;
end % for function inputIntermediateBearingHertz()




