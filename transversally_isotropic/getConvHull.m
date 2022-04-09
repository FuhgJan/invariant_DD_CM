function [invHull] = getConvHull(n_points, delta,a0)
% invHull: Establish points for convex hull of principal invariants as
%          outlined in Algorithm 1
% 
% arguments: (input)
%  n_points: Invariant set we would like to find a corresponding C to
%
%  delta: Delta of the deformation gradient as defined by eq.(45) in the paper
%
%  a0: Preferred direction of trans. isotropic material
%
% arguments: (output)
%  invHull: n_points points establishing the convex hull 

% Undeformed configuration of F in vector-form
mb = [1;0;0;0;1;0;0;0;1];

% Upper and lower bounds to sample from
lb = mb-delta;
ub = mb+delta;

% Generate n_points samples in the deformation gradient domain using LHS
F_inp =  lhs_scaled(n_points,lb,ub);
for j=1:n_points
    F = [F_inp(j,1), F_inp(j,2), F_inp(j,3);
        F_inp(j,4), F_inp(j,5), F_inp(j,6);
        F_inp(j,7), F_inp(j,8), F_inp(j,9)];
    
    C = F'*F;
    invHull(j,:) = getInvariants(C,a0);

end

end

