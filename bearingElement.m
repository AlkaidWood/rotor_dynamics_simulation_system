%% bearingElement
% generate the stiffness, damping matrix of a bearing element without mass
%% Syntax
% [Ke, Ce] = bearingElement(ANBearing)
%% Description
% ANBearing is a struct saving the physical parameters of a bearing element
% with fields: dofOfEachNodes, stiffness, damping
%
% Ke, Ce are stiffness, damping matrix of a bearing element. 
% (n*n, n is the number of dofs on this bearing element)
%% Symbols
% k: stiffness of bearing
%
% c: damping of bearing


function [Ke, Ce] = bearingElement(ANBearing)

% constants
k = ANBearing.stiffness;
c = ANBearing.damping;
dof = ANBearing.dofOnShaftNode;

%%

% generate stiffness matrix of bearing element
Ke = [ k, 0;
       0, k ];
Ke = blkdiag( Ke,zeros(dof - length(Ke)) ); % expand stiffness matrix

%%

% generate damping matrix of bearing element
Ce = [ c, 0;
       0, c ];
Ce = blkdiag( Ce,zeros(dof - length(Ce)) ); % expand damping matrix

end