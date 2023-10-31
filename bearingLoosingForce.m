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
 
 
if qn(66) >= 0 && qn(66) <= 0.0005
    K = Matrix.stiffnessLoosing;
    C = Matrix.dampingLoosing;
else
    K = Matrix.stiffness;
    C = Matrix.damping;
end
 
end
 
