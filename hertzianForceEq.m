%% herzianForceEq
% calculate the Herzian contact force for the roller bearing
%% Syntax
% f = herzianForceEq(tn, x, y, omegai, omegao, nb, ri, ro, delta0, kHertz, n)
%% Description
% tn: is the n-th time (s)
% 
% x,y: are the displacement at the n-th time in the x,y direction
% 
% omegai, omegao: are the rotation velosity of the inner and outer race of
% the roller bearing
%
% nb: is the number of the roller (ball or cylinder) in the bearing
%
% ri,ro: are radius of the inner and outer race of the bearing
%
% delta0: is the initial gap between roller and the inner race
%
% kHertz: is Herzian contact stiffness calculated by contact mechanics
% theory
%
% n: is the coefficient respect to the contact objects
%
% f = [fx; fy]: is a 2*1 vector representing the Herzian force in the x,y
% dirctions

function f = hertzianForceEq(tn, x, y, omegai, omegao, nb, ri, ro, delta0, kHertz, n)

omegac = (omegao*ro + omegai*ri) / (ro + ri); % rotation velosity of the cage
c1 = 2*pi/nb;
c2 = omegac*tn;
c3 = 0;
c4 = 0;
for ik = 1:1:nb
    thetak = c1 * (ik-1) + c2; % the angular position of the k-th roller in the bearing
    c5 = cos(thetak);
    c6 = sin(thetak);
    deltak = x*c5 + y*c6 - delta0; % the relative contact deformation between the k-th roller and the inner race of the bearing
    if deltak>0
        c7 = deltak^n;
        c3 = c3 + c7 * c5;
        c4 = c4 + c7 * c6;
    end % end if
end
f = kHertz * [c3; c4];

end