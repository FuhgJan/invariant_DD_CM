function Invariants = getInvariants(C,a0)
% getInvariants: Get the principal and pseudo invariants for trans. 
% isotropic materials of a right Cauchy-Green tensor and a preferred direction
% 
% arguments: (input)
%  C: Right Cauchy-Green tensor
%
% arguments: (output)
%  Invariants: I1,I2,I3,I4,I5
I1 = trace(C);
I2 = 0.5*(trace(C)^2 - trace(C*C));
I3 = det(C);

a0_kron_a0 = a0*a0';
I4 = trace(a0_kron_a0*C);
I5 = trace(a0_kron_a0*C'*C);

Invariants = [I1, I2, I3, I4, I5];

end

