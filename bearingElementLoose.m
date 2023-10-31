%% bearingElementLoose
% generate the stiffness, damping matrix of a bearing element with loosing
% fault
%% Syntax
% [Ke, Ce] = bearingElementLoose(AMBearing)
%% Description
% AMBearing is a struct saving the physical parameters of a bearing element
% with fields: dofOfEachNodes, stiffness, damping, mass, dofOnShaftNode,
% loosingStiffness, loosingdamping
%
% Ke, Ce are stiffness, damping cell of a bearing element. 
% (n*n, n is the number of dofs on this bearing element)
%% Symbols
% k: stiffness of bearing
%
% c: damping of bearing
%
% kl: stiffness of a loosing bearing
%
% cl: damping of a loosing bearing


function [Ke, Ce] = bearingElementLoose(AMBearing)

% constants
k = AMBearing.stiffness;
c = AMBearing.damping;
kl = AMBearing.loosingStiffness;
cl = AMBearing.loosingDamping;
dofShaft = AMBearing.dofOnShaftNode;
dofBearing = AMBearing.dofOfEachNodes;

%%

% stiffness matrix
Kin1 = [ k, 0;...
         0, k ];
Kin2 = [ 2*k, 0;...
         0,   k + kl ];
 
K11 = zeros(dofShaft);      K12 = zeros(dofShaft, dofBearing);
K21 = K12';                 K22 = zeros(dofBearing);

K11 = addElementIn(K11,Kin1,[1,1]);  K12 = addElementIn(K12,-Kin1,[1,1]);
K21 = addElementIn(K21,-Kin1,[1,1]); K22 = addElementIn(K22,Kin2,[1,1]);

Ke = {K11, K12;...
      K21, K22 };

%%

% damping matrix
Cin1 = [ c, 0;...
         0, cl ];
Cin2 = [ 2*c, 0;...
         0,   c + cl ];
 
C11 = zeros(dofShaft);      C12 = zeros(dofShaft, dofBearing);
C21 = C12';                 C22 = zeros(dofBearing);

C11 = addElementIn(C11,Cin1,[1,1]);  C12 = addElementIn(C12,-Cin1,[1,1]);
C21 = addElementIn(C21,-Cin1,[1,1]); C22 = addElementIn(C22,Cin2,[1,1]);

Ce = {C11, C12;...
      C21, C22 };

end  