%% hertzianForce
% calculate the Herizan contact force in the bearing
%% Syntax
% f = herzianForce(qn, tn, domega)
%% Description
% hertzDof: is an index n*2 matrix. The first column saves the dof No. of
% the Inner shaft, the second column saves the dof No. of the Outer shaft.
% If this functioin is used to calculte the Hertzian force of a Bearing
% connecting the ground, the second column would expect a dofNum+1. The
% i-th bearing info should be saved in the i-th row.
%
% omegaiNo, omegaoNo: are index vector denoting the inner/outer shaft No.
%
% n: is a vector saving the exponential power of the bearing. 
%
% kHertz: is a vector saving the Hertzian contact stiffness.
% 
% delta0: is a vector saving the initial gap of bearings
%
% ro, ri: is a vector saving the outer/inner radius of bearings
%
% nb: is a vector saving the nubmer of rollers in the bearing
%
% domega: is a vector saving the rotational speed of shafts
%
% tn: is a scalar denoting the time (s)
% 
% qn: is a vector saving the displacement of dofs
% 
% fHertz: is a dofNum*1 vector denoting the Hertzian force
function fHertz = hertzianForce(qn, tn, domega)
 
nb = [8   8   8  10   8   8];
ri = [0.0288      0.0288      0.0288      0.0755      0.0288      0.0288];
ro = [0.047       0.047       0.047        0.11       0.047       0.047];
delta0 = [1e-05       1e-05       1e-05     1.4e-05       1e-05       1e-05];
kHertz = [10800000000  10800000000  10800000000  14900000000  10800000000  10800000000];
n = [1.5         1.5         1.5         1.5         1.5         1.5];
omegaiNo = [1  1  1  2  1  1];
omegaoNo = [3  3  3  3  2  2];
hertzDof = [1   5  33  69  13  21;...
            107   77   81   91   95   65]';
 
domega = [domega, 0]; % add 1 row for adapting the bearing connecting the ground
qn = [qn; 0; 0]; % add 2 rows for adapting the bearing connecting the ground
 
fHertz2 = zeros(108,1); % add 2 rows for adapting the bearing connecting the ground
for iHertz = 1:1:6
    x = qn(hertzDof(iHertz,1)) - qn(hertzDof(iHertz,2)); % displacement of the Inner shaft - that of the Outer shaft
    y = qn(hertzDof(iHertz,1)+1) - qn(hertzDof(iHertz,2)+1);
    f = hertzianForceEq(tn, x, y, domega(omegaiNo(iHertz)), domega(omegaoNo(iHertz)), nb(iHertz), ri(iHertz), ro(iHertz), delta0(iHertz), kHertz(iHertz), n(iHertz));
    fHertz2(hertzDof(iHertz,1):hertzDof(iHertz,1)+1) = -f; % the relative displacement is (x_Inner - x_Outer), so the force on the Inner shaft should add a minus "-".
    fHertz2(hertzDof(iHertz,2):hertzDof(iHertz,2)+1) = f; % the force on the outer shaft/ground
end
fHertz = fHertz2(1:end-2);
 
end
 
