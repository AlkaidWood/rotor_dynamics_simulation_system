%% inputRubImpact
% Input parameters about rub impact fault
%% Syntax
%  OutputParameter = inputRubImpact(InputParameter)
%% Description
% InputParameter is a struct data saving all initial parameters except 
% parameters about rub impact fault.
% 
% OutputParameter is a strct data recording all InputParameter and the data
% about rub impact fault.
%
% All of parameters about rub impact should be typed into this .m file.


%%
function OutputParameter = inputRubImpact(InputParameter)

RubImpact.amount                = 1;
RubImpact.inShaftNo             = 1;
RubImpact.positionOnShaftDistance = 0.733; % m
RubImpact.interval              = 6.55e-3; % the gap between rotator and stator (m)
RubImpact.stiffness             = 1e6;     % N*m
RubImpact.coefficientOfFriction = 0.2;

checkInputData(RubImpact);

%%

OutputParameter = InputParameter;
OutputParameter.RubImpact = RubImpact;
OutputParameter.ComponentSwitch.hasRubImpact = true;
end