%% inputEssentialParameter
% Input parameters about shaft(s), disk(s), linear bearings and running
% status
%% Syntax
%  inputEssentialParameter()
%% Description
% There is no input parameter.
% 
% All essential initial parameters should be typed into this .m file.
% 
% All vectors in this .m file should be in column.


%%
function InitialParameter = inputEssentialParameter()

% typing the parameter about shaft
Shaft.amount            = 2;
Shaft.totalLength       = [0.85; 0.382]; % all vectors in column (m)
Shaft.dofOfEachNodes    = 4 * ones(Shaft.amount,1);
Shaft.outerRadius       = [0.01; 0.0325]; % m
Shaft.innerRadius       = [0; 0.02]; % m
Shaft.density           = 7890 * ones(Shaft.amount,1); % kg/m^3
Shaft.elasticModulus    = 209e9 * ones(Shaft.amount,1); % Pa
Shaft.poissonRatio      = 0.269 * ones(Shaft.amount,1);
checkInputData(Shaft)
Shaft.rayleighDamping   = [10, 0]; % [alpha, beta] CShaft = alpha*(MShaft+MDisk) + beta*KShaft

%%

% typing the parameter about running status
Status.ratio            = [-1.3]; % [v-shaft2/v-shaft1; v-shaft3/v-shaft1]
Status.vmax             = 500; % rad/s, the maximum rotational speed for shaft 1
Status.acceleration     = 100; % rad/s^2, acceleration of shaft 1
Status.duration         = 100; % s, the duration of shaft 1 in vmax
Status.isDeceleration   = false; % boolean, add a deceleration in status
Status.vmin             = 0; % s, the minimum speed afterdeceleration

% check input
if  length(Status.ratio) >= Shaft.amount
    error('too much input parameter in Status.ratio')
end

%%

% typing the parameter about disk
Disk.amount             = 4;
Disk.inShaftNo          = [1*ones(1,2), 2*ones(1,2)]'; % disks in the i-th shaft
Disk.dofOfEachNodes     = 4 * ones(Disk.amount,1);
Disk.radius             = [0.1215*ones(1,2), 0.124*ones(1,2)]'; % m
Disk.thickness          = [0.015*ones(1,4)]'; % m
Disk.positionOnShaftDistance = 1e-3 * [230.5, 664.5, 120.5, 225.5]'; %from left end (m)
Disk.density            = [7890*ones(1,4)]'; % kg/m^3
Disk.eccentricity       = [1e-5*ones(1,4)]'; % m

% check input
checkInputData(Disk)

for iDisk = 1:1:Disk.amount
   if Shaft.dofOfEachNodes(Disk.inShaftNo(iDisk)) ~= Disk.dofOfEachNodes(iDisk)
      error(['the dof of each disk should equal to the dof of the shaft'...
            ,' this disk locating']); 
   end
end

%%

% typing the parameter about linear bearing
% If you choose to input the bearing parameter here, you should not use the
% inputBearingHertz()
% model: shaft--k1c1--mass--k2c2--basement
% Bearing.amount          = 3;
% Bearing.inShaftNo       = [1; 1; 2];
% Bearing.dofOfEachNodes  = [2; 2; 2]; % if mass=0, dof must be 0 
% Bearing.positionOnShaftDistance = 1e-3 * [176.5; 718.5; 343.5];
% % stiffness = [bearing1_k1, bearing1_k2; bearing2_k1, bearing2_k2]
% Bearing.stiffness       = [1e8, 1e8; 1e8, 1e8; 1e8, 1e8]; % N*m
% % damping = [bearing1_c1, bearing1_c2; bearing2_c1, bearing2_c2]
% Bearing.damping         = [300, 300; 300, 300; 300, 300]; % N*s/m
% Bearing.mass            = [3; 3; 3]; % kg

Bearing.amount          = 0;
Bearing.inShaftNo       = [];
Bearing.dofOfEachNodes  = []; % if mass=0, dof must be 0 
Bearing.positionOnShaftDistance = [];
Bearing.stiffness       = []; % N*m
Bearing.damping         = []; % N*s/m
Bearing.mass            = []; % kg

checkInputData(Bearing)

%%

% ComponentSwitch will be changed by corresponding input..() function
ComponentSwitch.hasIntermediateBearing = false;
ComponentSwitch.hasLoosingBearing = false;
ComponentSwitch.hasRubImpact = false;
ComponentSwitch.hasCouplingMisalignment = false;
% ComponentSwitch.BearingH = false;
% ComponentSwitch.IntermediateBeaingH = false;

%%

% Output initialParameter without optional parameter
InitialParameter.Status          = Status;
InitialParameter.Shaft           = Shaft;
InitialParameter.Disk            = Disk;
InitialParameter.Bearing         = Bearing;
InitialParameter.ComponentSwitch = ComponentSwitch;


end