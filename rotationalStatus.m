%% rotationalStatus
% calculate the phase, speed, acceleratation for a typical rotation
%% Syntax
% status = rotationalStatus(t, vmax, duration)
%
% status = rotationalStatus(t, vmax, duration, acceleration)
%
% status = rotationalStatus(t, vmax, duration, acceleration, isDeceleration)
%
% status = rotationalStatus(t, vmax, duration, acceleration, isDeceleration, vmin)
%% Description
% t: is a row vector (unit, s)
%
% vmax: is maximum speed (rad/s)
%
% duration: is the time machinery running with maximum speed (s)
%
% acceleration: is an absolute value in running stauts (ac-/de- celerate, 
% default = 0)
%
% isDeceleration: is a boolean value control the deceleration in running
% status (default = false)
% 
% vmin: is minmum speed after decelerating (rad/s, default = 0)
%
% status = [omega; domega, ddomega] is a 3*n mareix. 


function status = rotationalStatus(t, vmax, duration, acceleration, isDeceleration, vmin)

% default value
if nargin < 6
    vmin = 0;
end
    
if nargin < 5
    isDeceleration = false;
end

if nargin < 4
    acceleration = 0;
end

%%

% check input
if t(1)<0 || duration < 0 
    error('time and duration must larger than 0 or equal 0');
end

% if  vmax <= 0
%     error('vmax must larger than 0');
% end

if abs(vmin) > abs(vmax)
    error('vmin should smaller than vmax')
end

%%

% constants
dataNum = length(t);
a = acceleration;
omega   = zeros(1,dataNum);
domega  = zeros(1,dataNum);
ddomega = zeros(1,dataNum);

%%

% acceleration = 0 
if (a == 0) && (isDeceleration == false)
    ddomega = zeros(1,dataNum);
    domega = vmax .* ones(1,dataNum);
    omega = vmax .* t;
elseif (acceleration == 0) && (isDeceleration == true)
    error('when acceleration equal 0, isDeceleration must be false')
end

%%

% acceleration > 0
if (abs(a) > 0) && (isDeceleration == false)
    t1 = abs(vmax/a); % the first key point (a>0 -> a=0)
    phase1 = 0.5 * a * t1^2; % the phase at t1
    
    ddomega = a                     .* (t <= t1)...
              + 0                   .* (t > t1);
          
    domega  = a .* t                .* (t <= t1)...
              + vmax                .* (t > t1);
          
    omega   = 0.5 * a .* t.^2                .* (t <= t1)...
              + ( phase1 + vmax .*(t-t1) )   .* (t > t1);
    
elseif (abs(a) > 0) && (isDeceleration == true)
    t1 = abs(vmax/a); % the first key point (a>0 -> a=0)
    phase1 = 0.5 * a * t1^2; % the phase at t1
    t2 = t1 + duration; % the end of uniform motion (a=0 -> a<0)
    phase2 = phase1 + vmax * (t2 - t1);
    t3 = t2 + (abs(vmax)-abs(vmin))/abs(a); % (a<0 -> a=0)
    phase3 = phase2 + vmax*(t3-t2) - 0.5*a*(t3-t2)^2;
    
    ddomega = a                     .* (t <= t1)...
              + 0                   .* ((t>t1) & (t<=t2))...
              - a                   .* ((t>t2) & (t<=t3))...
              + 0                   .* (t > t3);
          
    domega  = a .* t                .* (t <= t1)...
              + vmax                .* ((t>t1) & (t<=t2))...      
              + (vmax - a.*(t-t2))  .* ((t>t2) & (t<=t3))...
              + vmin                .* (t > t3);
          
    omega   = 0.5 * a .* t.^2                .* (t <= t1)...
              + (phase1 + vmax.*(t-t1))      .* ((t>t1) & (t<=t2))...
              + (phase2 + vmax.*(t-t2) - 0.5*a.*(t-t2).^2)  .* ((t>t2) & (t<=t3))...
              + (phase3 + vmin.*(t-t3))      .* (t > t3);
     
end % end if 

%%

% output
status = [omega; 
          domega; 
          ddomega];

end


