%% bearingElementMass
% generate the mass, stiffness, damping matrix of a bearing element with mass
%% Syntax
% [Me, Ke, Ce] = bearingElementMass(AMBearing)
%% Description
% AMBearing is a struct saving the physical parameters of a bearing element
% with fields: dofOfEachNodes, stiffness, damping, mass, dofOnShaftNode
%
% Me, Ke, Ce are mass, stiffness, damping cell of a bearing element. 
% (n*n, n is the number of dofs on this bearing element)
%% Symbols
% m: stiffness of bearing
% 
% k: stiffness of bearing
%
% c: damping of bearing


function [Me, Ke, Ce] = bearingElementMass(AMBearing)

% constants
k = AMBearing.stiffness;
c = AMBearing.damping;
m = AMBearing.mass;
dofShaft = AMBearing.dofOnShaftNode;
dofBearing = AMBearing.dofOfEachNodes;

%%

% stiffness matrix
Kin = [ k, 0;...
        0, k ];
 
K11 = zeros(dofShaft);      K12 = zeros(dofShaft, dofBearing);
K21 = K12';                 K22 = zeros(dofBearing);

K11 = addElementIn(K11,Kin,[1,1]);  K12 = addElementIn(K12,-Kin,[1,1]);
K21 = addElementIn(K21,-Kin,[1,1]); K22 = addElementIn(K22,2*Kin,[1,1]);

Ke = {K11, K12;...
      K21, K22 };

%%

% damping matrix
Cin = [ c, 0;...
        0, c ];
 
C11 = zeros(dofShaft);      C12 = zeros(dofShaft, dofBearing);
C21 = C12';                 C22 = zeros(dofBearing);

C11 = addElementIn(C11,Cin,[1,1]);  C12 = addElementIn(C12,-Cin,[1,1]);
C21 = addElementIn(C21,-Cin,[1,1]); C22 = addElementIn(C22,2*Cin,[1,1]);

Ce = {C11, C12;...
      C21, C22 };

%%

% mass matrix
Min = [ m, 0;...
        0, m ];
 
M11 = zeros(dofShaft);      M12 = zeros(dofShaft, dofBearing);
M21 = M12';                 M22 = zeros(dofBearing);

M22 = addElementIn(M22, Min, [1,1]);

Me = {M11, M12;...
      M21, M22 };
  
end  