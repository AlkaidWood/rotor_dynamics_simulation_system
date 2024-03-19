%% dynamicEquation
% saving the differential equation in this function
%% Syntax
% ddyn = dynamicEquation(tn, yn, dyn, Parameter)
%% Description
% tn: is the n-th time (s)
% 
% yn: is the displacement at the n-th time
% 
% dyn: is the first derivative at the n-th time
% 
% Parameter: is a struct saving all information about the model
% 
% ddyn: is the second derivative at the n-th time
 
 
function ddyn = dynamicEquation(tn, yn, dyn, Parameter)
 
 
% calculate phase, speed and acceleration
if tn <= 5
   ddomega = [100 -130];
   domega  = ddomega .* tn ;
   omega   = 0.5 * ddomega * tn^2;
else
   ddomega = zeros(1,2);
   domega  = [500 -650];
   omega   = [1250 -1625] + domega .* (tn-5);
end
 
 
% load matrix
M = Parameter.Matrix.mass;
G = Parameter.Matrix.gyroscopic;
N = Parameter.Matrix.matrixN;
Q = Parameter.Matrix.unblanceForce;
G(1:40, 1:40) = domega(1)*G(1:40, 1:40);
N(1:40, 1:40) = ddomega(1)*N(1:40, 1:40);
G(41:76, 41:76) = domega(2)*G(41:76, 41:76);
N(41:76, 41:76) = ddomega(2)*N(41:76, 41:76);
 
 
% calculate unblance force
diskInShaftNo = [1  1  2  2];
Q([9  29  53  61])   = [0.0054887   0.0054887   0.0057169   0.0057169] .* ( ddomega(diskInShaftNo) .* sin(omega(diskInShaftNo)) + domega(diskInShaftNo).^2 .* cos(omega(diskInShaftNo)));
Q([10  30  54  62]) = [0.0054887   0.0054887   0.0057169   0.0057169] .* (-ddomega(diskInShaftNo) .* cos(omega(diskInShaftNo)) + domega(diskInShaftNo).^2 .* sin(omega(diskInShaftNo)));
 
 
% calculate Hertzian force
fHertz = hertzianForce(yn, tn, domega);
 
 
% check the loosing bearing. If bearing loose, create the KLoose, CLoose.
[K, C] = bearingLoosingForce(yn, Parameter.Matrix);
 
 
% total force 
F = Q + fHertz;
 
 
% dynamic equation 
ddyn = M \ ( F -  (K - N)*yn - (C - G)*dyn );
 
end
 
