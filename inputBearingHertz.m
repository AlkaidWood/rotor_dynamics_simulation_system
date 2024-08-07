%% inputBearingHertz
% Input parameters about bearing considering Hertzian contact
%% Synatx
% Parameter = inputBearingHertz(InitialParameter)
%% Description 
% InitialParameter: is a struct generated by inputEssentialParameter()
%
% All parameters about bearing with Hertzian contact should be typed into
% this .m file.
%
% All vectors in this .m file should be in column
% 
% The contact stiffness could be calculated by main_contactStiffness.mlx


function Parameter = inputBearingHertz(InitialParameter)
% typing the parameters about bearing considering the Hertz contact between
% rollers and races
Bearing.amount          = 6;
Bearing.inShaftNo       = [1; 1; 2; 1; 1; 2];
Bearing.dofOfEachNodes  = [2,2,0; 2,2,0; 2,2,0; 0,0,0; 2,2,2; 0,0,0]; % if mass=0, dof must be 0 
Bearing.positionOnShaftDistance = 1e-3 * [176.5; 718.5; 343.5; 0; 850; 150];
Bearing.isHertzian      = [true; true; true; true; false; false]; % boolean
% M K C, elements in the same row: the MKC at the same position of the
% shaft; mass(1,1) -> mass(1,n):
% the mass near the shaft the bearing connecting ->
% the mass near the basement of bearing.
%
% If isHertizian and no mass, the corresponding k c will be added in
% global matrix normally; the model:
% shaft--Hertz+k1c1--basement;
% If is no Hertzian and with mass: there are n mass in a row, and n+1 k c 
% for a bearing; the model will be established as:
% shaft--k1c1--m1--k2c2--m2--k3c3--m3--k4c4- ...-mn--k(n+1)c(n+1)--basement;
% If is no Hertzian and no mass, the k c will be added in global 
% matrix normally; the model:
% shaft--k1c1--basement;
% If isHertzian and with mass, the hertzian force will be added at the mass
% in the first column (near the shaft); the model:
% shaft--Hertz+k1c1--m1--k2c2--m2--k3c3--m3--k4c4- ...-mn--k(n+1)c(n+1)--basement;

Bearing.stiffness       = [0,  2.5e7,  1e8, 0;...
                            0,  2.5e7,  1e8, 0;...
                            0,  2.5e7,  1e8, 0;...
                            1e6,    0,    0, 0;...
                            1e6,1.1e6,1.2e6, 1e6;...
                            1e6,    0,    0, 0]; % N*m
Bearing.damping         = [0, 150, 300, 0;...
                            0, 150, 300, 0;...
                            0, 150, 300, 0;...
                            20,  0,   0, 0;...
                            20,  10,  5, 10;...
                            20,  0,   0, 0]; % N*s/m
% the first n mass in each row must be non-zero, n is the number of mass
% of bearings
Bearing.mass            = [0.11, 3,    0;...
                            0.11, 3,    0;...
                            0.78, 2.5,  0;...
                            0,    0,    0;...
                            0.05, 0.05, 0.1;...
                            0,    0,    0]; % kg
Bearing.rollerNum = [8; 8; 10; 8; 0; 0];
Bearing.radiusInnerRace = [28.8e-3; 28.8e-3; 75.5e-3; 28.8e-3; 0; 0]; % m
Bearing.radiusOuterRace = [47e-3; 47e-3; 110e-3; 47e-3; 0; 0]; % m
Bearing.clearance = [10e-6; 10e-6; 14e-6; 10e-6; 0; 0]; % m
Bearing.contactStiffness = [1.08e10; 1.08e10; 1.49e10; 1.08e10; 0; 0]; % N*m^-3/2
Bearing.coefficient = [3/2; 3/2; 3/2; 3/2; 0; 0]; % =3/2 in a ball bearing; = 10/9 in a roller bearing

% check input data
checkInputData(Bearing)
% order the column in struct with shaft no and distance on the shaft
Bearing = sortRowsWithShaftDis(Bearing);

% Outpt initialParameter
Parameter = InitialParameter;
Parameter.Bearing = Bearing;
if sum(Bearing.isHertzian)~=0
    Parameter.ComponentSwitch.hasHertzianForce = true;
end % end if
end