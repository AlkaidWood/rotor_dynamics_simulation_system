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

LoosingBearing.amount = 1; % up to now, this program only support to simulate 1 loosing fault; amount = 1
% inbearingNo corresponding to the row index of the Parameter.Bearing.xxx after establishModel()
LoosingBearing.inBearingNo = [6]; 
LoosingBearing.interval = [1e-4]; % m
LoosingBearing.loosingStiffness = [1.5e7]; % N*m
LoosingBearing.loosingDamping = [100]; % N*s/m

% the model of bearing element:
% shaft--k1c1--m1--k2c2--m2--k3c3--m3--k4c4- ...-mn--k(n+1)c(n+1)--basement;
LoosingBearing.loosingPositionNo = [2]; % indicates the k,c no. which is loosing (e.g., loosing k1,c1 -> 1 or loosing k2,c2 -> 2)


checkInputData(LoosingBearing);

if InputParameter.Bearing.mass(LoosingBearing.inBearingNo) == 0
    error('the mass of Bearing which is loose must be lager than 0');
end

%%

OutputParameter = InputParameter;
OutputParameter.LoosingBearing = LoosingBearing;
OutputParameter.ComponentSwitch.hasLoosingBearing = true;

end