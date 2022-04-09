function Invariants = getInvariants(C)
% getInvariants: Get the principal invariants of a right Cauchy-Green tensor
% 
% arguments: (input)
%  C: Right Cauchy-Green tensor
%
% arguments: (output)
%  Invariants: I1,I2,I3

I1 = trace(C);
I2 = 0.5*(trace(C)^2 - trace(C*C));
I3 = det(C);

Invariants = [I1, I2, I3];

end

