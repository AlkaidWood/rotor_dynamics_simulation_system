%% misalignmentForce
% saving the equation of coupling misalignment force in this function
%% Syntax
% fMisalignment = misalignmentForce(qn, omega, domega)
%% Description
% omega, domega: are vector denoting the phase and speed of the shaft 
% 
% fMisalignment: is misalignment force (vector)
 
 
function fMisalignment = misalignmentForce(omega, domega)
 
 
fMisalignment = zeros(70,1);
for iMis = 1:1:1
    inShaftNo = [1];
    constants = [-0.12] .* domega(inShaftNo).^2 ;
    fMisalignment([17]) = constants .* sin(2*omega(inShaftNo));
    fMisalignment([17]) = constants .* cos(2*omega(inShaftNo));
end
 
end
 
