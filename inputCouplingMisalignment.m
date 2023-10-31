%% inputCouplingMisalignment
% Input parameters about coupling misalignment fault
%% Syntax
%  OutputParameter = inputCouplingMisalignment(InputParameter)
%% Description
% InputParameter is a struct data saving all initial parameters except
% parameters about coupling misalignment.
% 
% OutputParameter is a strct data recording all InputParameter and the data
% about coupling misalignment.
%
% All of parameters about coupling misalignmnet fault should be typed into this .m file.


%%
function OutputParameter = inputCouplingMisalignment(InputParameter)

CouplingMisalignment.amount = 1;
CouplingMisalignment.inShaftNo = 1;
CouplingMisalignment.positionOnShaftDistance = 0.733; % m
CouplingMisalignment.parallelOffset = 10e-3; % m
CouplingMisalignment.angleOffset = 0; % rad
CouplingMisalignment.distance = 50e-3; % m
CouplingMisalignment.mass = 10; % kg


checkInputData(CouplingMisalignment);

%%

OutputParameter = InputParameter;
OutputParameter.CouplingMisalignment = CouplingMisalignment;
OutputParameter.ComponentSwitch.hasCouplingMisalignment = true;
end