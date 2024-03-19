%% bearingLoosingForce
% saving the equation of loosing bearing force in this function
%% Syntax
% [K, C] = bearingLoosingForce(qn)
%% Description
% qn: is the displacement at the n-th time
% 
% Matrix: is a struct saving all matrix used in dynamic equation
% 
% K, C: are matrix after loosing
 
 
function [K, C] = bearingLoosingForce(qn, Matrix)
 
 
if qn(92) >= 0 && qn(92) <= 0.0001
    K = Matrix.stiffnessLoosing;
    C = Matrix.dampingLoosing;
else
    K = Matrix.stiffness;
    C = Matrix.damping;
end
 
end
 
