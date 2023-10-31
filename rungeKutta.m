%% rungeKutta
% calculate the value of second order differential equation after one step
%% Syntax
% [yn1,dyn1] = rungeKutta(fun,tn,yn,dyn,h) 
%
% [yn1,dyn1] = rungeKutta(fun,tn,yn,dyn) 
%% Description
% fun: is a function of second order dynamic equation. 
%
% $\ddot{x} = f\left( t,x,\dot{x}\right)$
%
% tn: is independent variable in the n-th step
%
% yn: is denpendent variable in the n-th step
%
% dyn: is the first derivative of the dependent variable in the n-th step
%
% h: h is a scalar denoting the step of this calcualtion (default, 1)
%
% yn1: is denpendent variable in the (n+1)-th step
%
% dyn1: is the first derivative of the dependent variable in the (n+1)-th step


function [yn1,dyn1] = rungeKutta(fun,tn,yn,dyn,h) 

%Check inputs
if nargin < 5
    h=1;
    if nargin < 4
         error('rungeKutta(): not enough inputs')
    end
end

%%
%input K11~K24
k11 = dyn;
k21 = fun(tn,yn,dyn);
k12 = dyn+h/2*k21;
k22 = fun(tn+h/2, yn+h/2*k11, dyn+h/2*k21);
k13 = dyn+h/2*k22;
k23 = fun(tn+h/2, yn+h/2*k12, dyn+h/2*k22);
k14 = dyn+h*k23;
k24 = fun(tn+h, yn+h*k13, dyn+h*k23);
yn1 = yn + h/6*(k11+ 2*k12+ 2*k13 +k14);
dyn1 = dyn + h/6*(k21 +2*k22 +2*k23 +k24);

end


