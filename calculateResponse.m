%% calculateResponse
% calculate the time history of the rotor system
%% Syntax
% [q, dq, t] = calculateResponse(Parameter, tSpan, samplingFrequency)
%% Description
% Parameter: is a struct saving the model data
%
% tSpan = [tStart, tEnd] is a 1*2 array saving the solving information
%
% samplingFrequency: is a integer equaling to 1/step
%
% isPlotStatus: is boolean to contol the plot of running status
%
% reduceInterval: denotes the re-sampling interval (scaler)
%
% q is time history of displacemnet (2D matrix, node * tNum)
%
% dq is time history of velocity (2D matrix, node * tNum)
% 
% t is time series (row)
%
% convergenceStr: is a str deliver the message of convergence


function [q, dq, t, convergenceStr] = calculateResponse(Parameter, tSpan, samplingFrequency, isPlotStatus, reduceInterval)

% check input
if nargin < 5
    reduceInterval = 1;
end

if nargin < 4
    isPlotStatus = true;
end

%%

% generate time series
tStart = tSpan(1);
tEnd = tSpan(2);
tNum = floor((tEnd - tStart) * samplingFrequency);
step = 1/samplingFrequency;
t = linspace(tStart,tEnd,tNum);

%%

% initial the response
dofNum   = Parameter.Mesh.dofNum;
q = zeros(dofNum, tNum); % for saving response
dq = zeros(dofNum, tNum);
yn = zeros(dofNum,1); % initial condition of differential euqtion
dyn = zeros(dofNum,1);
equation = @(tn,yn,dyn)dynamicEquation(tn,yn,dyn,Parameter);
convergenceStr = [];
% calculate response
for iT = 1:1:tNum
    [q(:,iT), dq(:,iT)] = rungeKutta(equation, t(iT), yn, dyn, step);
    yn = q(:,iT);
    dyn = dq(:,iT);
    if isnan(dyn)
        convergenceStr = ['t=', num2str(t(iT)), 's non-convergent'];
        break
    end
end

%%

%reduce data (re-sampling)
reduce_index = 1 : reduceInterval : size(q, 2); % index of the sampling data
if reduce_index(end) ~= size(q, 2)
   reduce_index = [reduce_index, size(q, 2)]; 
end
q  = q(:, reduce_index);% re_sampling
dq = dq(:,reduce_index);
t  = t(:,reduce_index);

%% 

% calculate the acceleration, rotational velocity and phase
if isPlotStatus
    % load constants
    Status   = Parameter.Status;
    shaftNum = Parameter.Shaft.amount;
    dofNum   = Parameter.Mesh.dofNum;
    vmax            = Status.vmax;
    duration        = Status.duration;
    acceleration    = Status.acceleration;
    isDeceleration  = Status.isDeceleration;
    vmin            = Status.vmin;
    ratio           = Status.ratio;
    ratio           = [1; ratio]; % the first shaft is basic


    % initial the status matrix
    omega 	= zeros(dofNum,length(t));
    domega  = omega;
    ddomega = omega;


    % calculate the status matrix
    status = cell(shaftNum,1);
    for iShaft = 1:1:shaftNum
        iVmax           = vmax * ratio(iShaft);
        iVmin           = vmin * ratio(iShaft);
        iAcceleration   = acceleration*ratio(iShaft);
        status{iShaft} = rotationalStatus(t, iVmax, duration, iAcceleration,...
                                          isDeceleration, iVmin);
    end


    % assemble the status matrix
    Node = Parameter.Mesh.Node;
    dofOnNodeNo = Parameter.Mesh.dofOnNodeNo;
    for iDof = 1:1:dofNum
        isBearing = Node( dofOnNodeNo(iDof) ).isBearing;
        if ~isBearing
            shaftNo = Node( dofOnNodeNo(iDof) ).onShaftNo;
            omega(iDof,:)   = status{shaftNo}(1,:);
            domega(iDof,:)  = status{shaftNo}(2,:);
            ddomega(iDof,:) = status{shaftNo}(3,:);
        end
    end


    % plot running status
    plotRunningStatus(t,status)
end


end