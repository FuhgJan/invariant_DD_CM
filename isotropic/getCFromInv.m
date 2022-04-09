function [C] = getCFromInv(I1,I2,I3)
% inhull: From a set of isotropic invariants deduce a possible right
% Cauchy Green tensor (eq. (46-49) in paper)
%
% arguments: (input)
%  I1,I2,I3: Principal invariant set we want to find a right Cauchy-Green
%            tensor for
%
% arguments: (output)
%  C: Possible right Cauchy-Green tensor


I = [I1,I2,I3];

% Eq. (46)
H = (1/9)*(I(1)^2 - 3*I(2));
G = (1/3)*I(1)*I(2) - I(3) - 2*(I(1)/3)^3;
if H==0 && G ==0
beta = pi/2;
else
beta = acos(-G/(2*(H)^(3/2)));
end

% Eq. (48)
lam31_sq = I(1)/3 + 2*(H)^(1/2)*cos(beta/3);
lam21_sq = I(1)/3 - 2*(H)^(1/2)*cos((pi+beta)/3);
lam11_sq = I(1)/3 - 2*(H)^(1/2)*cos((pi-beta)/3);

% Check if input invariance actually correspond to a physical state
test = (imag(lam31_sq)<1e-4) & (imag(lam21_sq)<1e-4) & (imag(lam11_sq)<1e-4);
if numel(find(test==0))>0
    error('Not real')
end

% Eq. (49)
C =  [real(lam11_sq),0,0;
    0,real(lam21_sq),0;
    0,0,real(lam31_sq)];

% Here we could check if what we found is correct
I1Test = trace(C);
I2Test = 0.5*(trace(C)^2- trace(C*C));
I3Test = det(C);

end

