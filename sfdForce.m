%% sfdForce
% saving the equation of the nonlinear force form sfd
%% Syntax
% fSfd = sfdForce(qn, dqn)
%% Description
% qn: is the displacement at the n-th time
% 
% dqn: is the speed at the n-th time 
%
% fSfd: is nonliear sfd force (vector)

function fSfd = sfdForce(qn, dqn)
fSfd = zeros(16,1);

positionOnShaftNode = 3;
radius = 53.3e-3; % m
mu = 0.22; % N*s/m
c = 0.4e-3; % clearance 0.35
l = 30e-3;


x  =  qn(positionOnShaftNode*4 - 3);
y  =  qn(positionOnShaftNode*4 - 2);
dx = dqn(positionOnShaftNode*4 - 3);
dy = dqn(positionOnShaftNode*4 - 2);

% % method 1-----------------------------------------------------------------
% constant = -mu*radius*l^3 / c^2;
% e = sqrt(x^2+y^2);
% r = e/c;
% dr = (x*dx + y*dy) / (c*e + eps);
% dpsi = (-1*y*dx + x*dy) / (x^2 + y^2 + eps);
% 
% 
% if r == 0
%     fx = 0;
%     fy = 0;
% else
%     theta1 = atan(-dr/(r*dpsi+eps));
%     theta2 = theta1 + pi;
%     
%     intEq11 = @(theta) sin(theta).*cos(theta)./(1+r.*cos(theta)).^3;
%     intEq02 = @(theta) cos(theta).^2./(1+r.*cos(theta)).^3;
%     intEq20 = @(theta) sin(theta).^2./(1+r.*cos(theta)).^3;
%     
%     int11 = integral(intEq11, theta1, theta2);
%     int02 = integral(intEq02, theta1, theta2);
%     int20 = integral(intEq20, theta1, theta2);
%     
%     fr = constant*(int11*dpsi*r + int02*dr);
%     ft = constant*(int20*dpsi*r + int11*dr);
%     
%     fx = x/e*fr - y/e*ft;
%     fy = y/e*fr + x/e*ft;
% end % end if
% %--------------------------------------------------------------------------

% method 2-----------------------------------------------------------------
constant = -mu*radius*l^3 ;
frac1 = 2*(c^2-x^2-y^2)^(3/2);
frac2 = 2*(c^2-x^2-y^2)^(5/2);
fx = constant*( (pi*dx)/frac1 +  (3*pi*x*(x*dx+y*dy))/frac2 );
fy = constant*( (pi*dy)/frac1 +  (3*pi*y*(x*dx+y*dy))/frac2 );

% -------------------------------------------------------------------------

fSfd(positionOnShaftNode*4 - 3) = fx;
fSfd(positionOnShaftNode*4 - 2) = fy;


end % end function