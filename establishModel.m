%% establishModel
% Generate the matrices and vectors in dynamic equation of general rotor system
%% Syntax
% Parameter = establishModel(InitialParameter)
% Parameter = establishModel(InitialParameter,NameValueArgs)
%% Description
% InitialParameter:
%                    Shaft: [1×1 struct]
%                     Disk: [1×1 struct]
%                  Bearing: [1×1 struct]
%          ComponentSwitch: [1×1 struct]
%      IntermediateBearing: [1×1 struct]
%                RubImpact: [1×1 struct]
%           LoosingBearing: [1×1 struct]
%     CouplingMisalignment: [1×1 struct]
%
% Parameter:
%                    Shaft: [1×1 struct]
%                     Disk: [1×1 struct]
%                  Bearing: [1×1 struct]
%          ComponentSwitch: [1×1 struct]
%      IntermediateBearing: [1×1 struct]
%                RubImpact: [1×1 struct]
%           LoosingBearing: [1×1 struct]
%     CouplingMisalignment: [1×1 struct]
%                     Mesh: [1×1 struct]
%                   Matrix: [1×1 struct]
%
% NameValueArgs.gridFineness = 'low' (default) or 'middle' or 'high' or
% cell data (manual grid)
%
% NameValueArgs.isPlotModel = true (default) or false
%
% NameValueArgs.isPlotMesh = true (default) or false



function Parameter = establishModel(InitialParameter,NameValueArgs)

arguments % name value pair
    InitialParameter 
    NameValueArgs.isPlotModel = true;
    NameValueArgs.isPlotMesh = true;
    NameValueArgs.gridFineness = 'low';
end

%%

% plot the schematic diagram of the rotor
if NameValueArgs.isPlotModel
    plotModel(InitialParameter);
end

%%

% mesh
Parameter = meshModel(InitialParameter,NameValueArgs.gridFineness); % choose: low, middle, high 
if NameValueArgs.isPlotMesh
	plotMesh(Parameter); % plot mesh result
end

%%

% generate FEM matrices of shaft
[MShaft, KShaft, GShaft, NShaft] = femShaft( Parameter.Shaft,...
                                             Parameter.Mesh.nodeDistance );

%%

% generate FEM matrices of disk 
[MDisk, GDisk, NDisk, QDisk] = femDisk( Parameter.Disk,...
                                        [Parameter.Mesh.Node.dof] );
                             

%%

% generate FEM matrices of bearing
if ~InitialParameter.ComponentSwitch.hasLoosingBearing
    [MBearing, KBearing, CBearing] = femBearing( Parameter.Bearing,...
                                               [Parameter.Mesh.Node.dof] );
else
    [MBearing, KBearing, CBearing, KLBearing, CLBearing] = femBearing( ...
                                                Parameter.Bearing,...
                                               [Parameter.Mesh.Node.dof],...
                                                Parameter.LoosingBearing);
end
                                         


%%

% generate FEM matrices of intermediate bearing
KInterBearing = zeros( length(MDisk) );
CInterBearing = zeros( length(MDisk) );
if Parameter.ComponentSwitch.hasIntermediateBearing
    [MInterBearing, KInterBearing, CInterBearing] = femInterBearing( ...
            Parameter.IntermediateBearing, [Parameter.Mesh.Node.dof] );
end


%%

% expand the matrices about shaft
shaftDofNum = length(MShaft);
diskDofNum = length(MDisk);
if diskDofNum > shaftDofNum
    MShaft = blkdiag( MShaft,zeros(diskDofNum - shaftDofNum) );
    KShaft = blkdiag( KShaft,zeros(diskDofNum - shaftDofNum) );
    GShaft = blkdiag( GShaft,zeros(diskDofNum - shaftDofNum) );
    NShaft = blkdiag( NShaft,zeros(diskDofNum - shaftDofNum) );
end
clear shaftDofNum diskDofNum;

%%

% rayleigh damping
rayleighCoeff = Parameter.Shaft.rayleighDamping;
CShaft = rayleighCoeff(1) * (MShaft+MDisk) + rayleighCoeff(2) * KShaft;


% assemble
M = MShaft + MDisk + MBearing + MInterBearing;
K = KShaft +         KBearing + KInterBearing;
G = GShaft + GDisk;
N = NShaft + NDisk;
C = CShaft +         CBearing + CInterBearing;
Q = QDisk;

if InitialParameter.ComponentSwitch.hasLoosingBearing
    KLoosing = KShaft + KLBearing + KInterBearing;
    CLoosing = CShaft + CLBearing + CInterBearing;
    KLoosing = sparse(KLoosing);
    CLoosing = sparse(CLoosing);
end


% sparse format
M = sparse(M);
K = sparse(K);
C = sparse(C);


% generate the struct Matrix
Matrix.mass = M;
Matrix.stiffness = K;
Matrix.gyroscopic = G;
Matrix.damping = C;
Matrix.matrixN = N;
Matrix.unblanceForce = Q;

if InitialParameter.ComponentSwitch.hasLoosingBearing
    Matrix.stiffnessLoosing = KLoosing;
    Matrix.dampingLoosing = CLoosing;
end

%%

% output

Parameter.Matrix = Matrix;

end