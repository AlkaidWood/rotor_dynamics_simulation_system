%% diskElement
% generate the mass, stiffness, gyroscopic matrix of a disk element
%% Syntax
% [Me, Ge, Ne] = diskElement(ADisk)
%% Description
% ADisk is a struct saving the physical parameters of a disk element with
% fields: dofOfEachNodes, radius, density, thickness, eccentricity
%
% Me, Ge, Ne are mass, gyroscopic and N matrix of a disk element. 
% (n*n, n is the number of all dofs on this shaft element)
%% Symbols
% m: mass of disk
%
% r: radius of disk
%
% rho: density 
%
% dof: degree of freedom of this disk element
%
% Id: rotational inertial about diameter
% 
% Ip: polar rotational inertial
%
% eDisk: eccentricity of the disk


function [Me, Ge, Ne] = diskElement(ADisk)

% check the input
fieldName = {'dofOfEachNodes', 'radius', 'density', 'thickness'};
hasFieldName = isfield(ADisk, fieldName);
if length(hasFieldName) ~= sum(hasFieldName)
    error('Incorrect field names for input struct')
end

%%

% calculate the constants
r = ADisk.radius;
thickness = ADisk.thickness;
rho = ADisk.density;
eDisk = ADisk.eccentricity;
m = (pi*r^2) * thickness * rho;
Id = 1/4 * m * r^2;
Ip = 1/2 * m * r^2;

%%

% mass matrix
MT = [ m, 0, 0, 0;...
       0, m, 0, 0;...
       0, 0, 0, 0;...
       0, 0, 0, 0 ];

MR = [ 0,  0,  0,  0;...
       0,  0,  0,  0;...
       0,  0, Id,  0;...
       0,  0,  0, Id ]; 
   
Me = MT + MR;

%%

% gyrosocpic matrix
Ge = [  0,   0,   0,   0;...
        0,   0,   0,   0;...
        0,   0,   0, -Ip;...
        0,   0,  Ip,   0 ]; 

   
%%

% Ne matrix
Ne = [  0,   0,   0,   0;...
        0,   0,   0,   0;...
        0,   0,   0,   0;...
        0,   0,  Ip,   0 ]; 
    

end