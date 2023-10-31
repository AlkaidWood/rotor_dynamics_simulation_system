%% rubImpactForce
% saving the equation of rub-impact force in this function
%% Syntax
% fRub = rubImpactForce(qn)
%% Description
% qn: is the displacement at the n-th time
% 
% fRub: is rub-impact force (vector)
 
 
function fRub = rubImpactForce(qn)
 
kRub   = [1000000];
mu     = [0.2];
delta  = [0.00655];
rubDof = [17];
 
 
fRub = zeros(70,1);
for iRub = 1:1:1
    e = sqrt(qn(rubDof(iRub))^2 + qn(rubDof(iRub)+1)^2);
    if e >= delta(iRub)
       fRub(rubDof(iRub))   = kRub(iRub)*(1-delta/e) * (qn(rubDof(iRub)) - mu*qn(rubDof(iRub)+1));
       fRub(rubDof(iRub)+1) = kRub(iRub)*(1-delta/e) * (mu*qn(rubDof(iRub)) + qn(rubDof(iRub)+1));
    end
end
 
end
 
