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
function OutputParameter = inputIntermediateBearing(InputParameter)

% typing the parameter about intermediate bearing
IntermediateBearing.amount          = 4;
% shaft no. connected by same bearing in row; different bearings in column
IntermediateBearing.betweenShaftNo  =  [2, 1;...
                                        1, 2;...
                                        1, 2;...
                                        1, 2]; % n*2
% dof
IntermediateBearing.dofOfEachNodes =  [2,0,0;...
                                       2,2,0;...
                                       2,2,2;...
                                       0,0,0];% if mass=0, dof must be 0 
% the same bearing in row; different bearings in column, n*2
IntermediateBearing.positionOnShaftDistance = 1e-3 *  [0,   258;...
                                                       400, 50;...
                                                       600, 100;...
                                                       500, 300];
IntermediateBearing.isHertzian      = [true;false;false;true]; % boolean
IntermediateBearing.isHertzianTop   = [false; false; false; true];
% M K C, elements in the same row: the MKC at the same position of the
% shaft; mass(1,1) -> mass(1,n): 
% the mass near the betweenShaftNo(:,1) -》the mass near the betweenShaftNo(:,2)
% If isHertizian and no mass, the corresponding k c will be added in
% global matrix normally; the model:
% shaft1--Hertz+k1c1--shaft2;
% If is no Hertzian and with mass: there are n mass in a row, and n+1 k c 
% for a bearing; the model will be established as:
% shaft1--k1c1--m1--k2c2--m2--k3c3--m3--k4c4- ...-mn--k(n+1)c(n+1)--shaft2;
% If is no Hertzian and no mass, the k c will be added in global 
% matrix normally; the model:
% shaft1--k1c1--shaft2;
% If isHertzian and with mass, the hertzian force will be added at the mass
% in the first column (near the shaft); the model:
% shaft1--Hertz+k1c1--m1--k2c2--m2--k3c3--m3--k4c4--mn--k(n+1)c(n+1)--shaft2; (isHertzianTop=true)
% shaft1--k1c1--m1--k2c2--m2--k3c3--m3--k4c4--mn--Hertz+k(n+1)c(n+1)--shaft2; (isHertzianTop=false)
IntermediateBearing.stiffness       =  [1e9, 0,   0,   0;...
                                        1e3, 1e3, 1e3, 0;...
                                        1e3, 1e3, 1e3, 1e3;...
                                        1e3, 0,   0,   0]; % N/m, in column, n*1
IntermediateBearing.damping         =  [300, 0,   0,   0;...
                                        50,  100, 150, 0;...
                                        50,  100, 150, 300;...
                                        300, 0,   0,   0]; % N/s^2, in column, n*1
IntermediateBearing.mass            =  [3, 0, 0;...
                                        3, 4, 0;...
                                        3, 4, 5;...
                                        0, 0, 0]; % kg
% if there is no Hertizan contact force, please set n*1 zero vector for following parameters, where n is the number of the intermediate bearing                                 
IntermediateBearing.rollerNum        = [8;0;0;8];
IntermediateBearing.radiusInnerRace = [28.8e-3;0;0;28.8e-3]; % m
IntermediateBearing.radiusOuterRace = [47e-3;0;0;47e-3]; % m
IntermediateBearing.clearance = [10e-6;0;0;10e-6]; % m
IntermediateBearing.contactStiffness = [1.08e10;0;0;1.08e10]; % N*m^-3/2
IntermediateBearing.coefficient = [3/2;0;0;3/2]; % =3/2 in a ball bearing; = 10/9 in a roller bearing

%%

%check the input data
if size(IntermediateBearing.betweenShaftNo, 2) ~= 2
    error('IntermediateBearingH.betweenShaftNo must be a 2 column matrix')
end
               
if size(IntermediateBearing.positionOnShaftDistance, 2) ~= 2
    error('IntermediateBearingH.positionOnShaftDistance must be a 2 column matrix')
end

checkInputData(IntermediateBearing);

%%

% sort the betweenShaftNo and PositionOnShaftDistance (sort column)
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
       IntermediateBearing.betweenShaftNo = shaftNo;
       IntermediateBearing.positionOnShaftDistance = position;
       % isHertzianTop, stiffness, damping, mass, dofOfEachNodes should be 
       % updated, if iBearing has mass.
       if sum(IntermediateBearing.mass(iBearing,:))
            % adjust isHertzianTop
            IntermediateBearing.isHertzianTop(iBearing) = ~IntermediateBearing.isHertzianTop(iBearing);
            % adjust mass and dofOfEachNodes
            colNum = size(IntermediateBearing.mass(iBearing,:),2); 
            massHere = IntermediateBearing.mass(iBearing,:);
            dofHere = IntermediateBearing.dofOfEachNodes(iBearing,:);
            massHereNum = sum(massHere~=0);
            massHere = massHere(1:massHereNum);
            dofHere = dofHere(1:massHereNum);
            massInv = flip(massHere);
            dofInv = flip(dofHere);
            IntermediateBearing.mass(iBearing,:) = [massInv, zeros(1,colNum-massHereNum)];
            IntermediateBearing.dofOfEachNodes(iBearing,:) = [dofInv, zeros(1,colNum-massHereNum)];
            % adjust stiffness and damping
            colNum = size(IntermediateBearing.stiffness(iBearing,:),2);
            kHere = IntermediateBearing.stiffness(iBearing,:);
            cHere = IntermediateBearing.damping(iBearing,:);
            kHereNum = massHereNum + 1;
            kHere = kHere(1:kHereNum);
            cHere = cHere(1:kHereNum);
            kInv = flip(kHere);
            cInv = flip(cHere);
            IntermediateBearing.stiffness(iBearing,:) = [kInv, zeros(1,colNum-kHereNum)];
            IntermediateBearing.damping(iBearing,:) = [cInv, zeros(1,colNum-kHereNum)];
       end
    end
end
% sort columns in struct
IntermediateBearing = sortRowsWithShaftDis(IntermediateBearing);

%%

OutputParameter = InputParameter;
OutputParameter.IntermediateBearing = IntermediateBearing;
OutputParameter.ComponentSwitch.hasIntermediateBearing = true;
if sum(IntermediateBearing.isHertzian)~=0
    OutputParameter.ComponentSwitch.hasHertzianForce = true;
end % end if
end % for function inputIntermediateBearingHertz()




