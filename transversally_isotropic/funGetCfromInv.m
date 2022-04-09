function [C] = funGetCfromInv(I,delta, a0)
% funGetCfromInv: Establish and solve the optimization problem of eq. (55) to find right
%                 Cauchy Green tensor to a specific set of transversally
%                 isotropic invariants
%
% arguments: (input)
%  I: Vector of set of invariants (I1,I2,I3,I4,I5)
%  
%  delta: \delta of the deformation gradient as defined by eq.(45) in the paper
%
%  a0: Preferred direction of trans. isotropic material
%
% arguments: (output)
%  C: Possible right Cauchy-Green tensor


% Definition of optimization procedure
fun = @(x)numSolveInv(I,x,a0);
options = optimoptions('fsolve','Display','none','StepTolerance', 1e-8,'FunctionTolerance', 1e-8);

% Bounds obtained from [F_{0}-delta, F_{0}+delta] where F_{0} is the 
% deformation gradient of the undeformed configuration 
lb = [(1-delta)^2;          
    2*(1+delta)*(-delta);
    2*(1+delta)*(-delta); 
    (1-delta)^2;
    2*(1+delta)*(-delta);
    (1-delta)^2];

ub = [(1+delta)^2+2*delta^2;
    2*(1+delta)*delta + delta^2;
    2*(1+delta)*delta + delta^2; 
    (1+delta)^2+2*delta^2;
    2*(1+delta)*delta + delta^2;
    (1+delta)^2+2*delta^2];

% Initial points
x0 = (lb+ub)/2;
% Solve optimization problem
[x,f] =fsolve(fun,x0, options);

% Maximum deviation from 0
disp(max(f));

% Possible C for the given invariant set I
C = [x(1) x(2) x(3);
    x(2) x(4) x(5);
    x(3) x(5) x(6)];

% Check if invariants of found C correspond to I
I_test = getInvariants(C,a0);

end

