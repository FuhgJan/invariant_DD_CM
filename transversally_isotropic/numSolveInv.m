function F = numSolveInv(I,x, a0)
% numSolveInv: Definition of optimization problem of eq. (55)
% 
% arguments: (input)
%  I: Invariant set we would like to find a corresponding C to
%
%  x: Values to be optimized (values of C)
%
%  a0: Preferred direction of trans. isotropic material
%
% arguments: (output)
%  F: Cost function
a0_kron_a0 = a0*a0';

C = [x(1) x(2) x(3);
     x(2) x(4) x(5);
     x(3) x(5) x(6)];
 
 
F(1) = trace(C) - I(1);

F(2) = 0.5*(trace(C)^2 - trace(C*C)) - I(2);
F(3) = det(C) - I(3);

F(4) = trace(a0_kron_a0*C)- I(4);
F(5)=  trace(a0_kron_a0*C*C) - I(5);

end
