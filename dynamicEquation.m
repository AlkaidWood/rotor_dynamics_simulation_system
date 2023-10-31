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
K = Parameter.Matrix.stiffness;
C = Parameter.Matrix.damping;
G(1:28, 1:28) = domega(1)*G(1:28, 1:28);
N(1:28, 1:28) = ddomega(1)*N(1:28, 1:28);
G(29:48, 29:48) = domega(2)*G(29:48, 29:48);
N(29:48, 29:48) = ddomega(2)*N(29:48, 29:48);
 
 
% calculate unblance force
diskInShaftNo = [1  1  2  2];
Q([9  17  33  37])   = [5.4887e-05  5.4887e-05  5.7169e-05  5.7169e-05] .* ( ddomega(diskInShaftNo) .* sin(omega(diskInShaftNo)) + domega(diskInShaftNo).^2 .* cos(omega(diskInShaftNo)));
Q([10  18  34  38]) = [5.4887e-05  5.4887e-05  5.7169e-05  5.7169e-05] .* (-ddomega(diskInShaftNo) .* cos(omega(diskInShaftNo)) + domega(diskInShaftNo).^2 .* sin(omega(diskInShaftNo)));
 
 
% total force 
F = Q;
 
 
% dynamic equation 
ddyn = M \ ( F -  (K - N)*yn - (C - G)*dyn );
 
end
 
