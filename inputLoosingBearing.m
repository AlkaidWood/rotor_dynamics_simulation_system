%% inputLoosingBearing
% Input parameters about basedment looseness fault 
%% Syntax
%  OutputParameter = inputLoosingBearing(InputParameter)
%% Description
% InputParameter is a struct data saving all initial parameters except
% parameters about basedment looseness fault.
% 
% OutputParameter is a strct data recording all InputParameter and the data
% about basedment looseness fault.
%
% All of parameters about basedment looseness fault should be typed into this .m file.


%%
function OutputParameter = inputLoosingBearing(InputParameter)

LoosingBearing.amount = 1; % up to now, this program support amount = 1;
LoosingBearing.inBearingNo = 1; 
LoosingBearing.interval = 5e-4; % m
LoosingBearing.loosingStiffness = 0.1e7; % N*m
LoosingBearing.loosingDamping =  400; % N*s/m


checkInputData(LoosingBearing);

if InputParameter.Bearing.mass(LoosingBearing.inBearingNo) == 0
    error('the mass of Bearing which is loose must be lager than 0');
end

%%

OutputParameter = InputParameter;
OutputParameter.LoosingBearing = LoosingBearing;
OutputParameter.ComponentSwitch.hasLoosingBearing = true;

end