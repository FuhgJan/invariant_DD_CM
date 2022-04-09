function [genPoints] = FunSimAnnealing45(genPoints,a0)
% FunSimAnnealing45: Use of simulated annealing to generate space-filling points in
%         isotropic invariant space (Algorithm 2 of paper)
%
% arguments: (input)
%  genPoints: Samples with evenly spaced primary invariant space
%
%  a0: Preferred direction of trans. isotropic material
%
% arguments: (output)
%  genPoints: Sampled points through the algorithm



%% Check physicality of generated points in principal invariant space and
% make sure that I4,I5 are a physically attainable set (eq. (46-49))
H = (1/9)*(genPoints(:,1).^2 - 3*genPoints(:,2));
G = (1/3)*genPoints(:,1).*genPoints(:,2) - genPoints(:,3) - 2*(genPoints(:,1)./3).^3;

beta = acos(-G./(2*(H).^(3/2)));
if H(1) ==0 && G(1)==0
    beta(1)=pi/2;
end
lam31_sq = genPoints(:,1)./3 + 2*(H).^(1/2).*cos(beta./3);
lam21_sq = genPoints(:,1)./3 - 2*(H).^(1/2).*cos((pi+beta)./3);
lam11_sq = genPoints(:,1)./3 - 2*(H).^(1/2).*cos((pi-beta)./3);

% Check numerical problems 
if numel(find(imag(lam31_sq)>1e-4)) ==0
    lam31_sq = real(lam31_sq);
end

if numel(find(imag(lam21_sq)>1e-4)) ==0
    lam21_sq = real(lam21_sq);
end

if numel(find(imag(lam11_sq)>1e-4)) ==0
    lam11_sq = real(lam11_sq);
end

C_mat = zeros(3,3,numel(lam11_sq));
genPoints_new = zeros(size(genPoints));
for i=1:numel(lam11_sq)
    C_mat(:,:,i) =  [lam11_sq(i),0,0;
        0, lam21_sq(i),0;
        0,0,lam31_sq(i)];
    
    genPoints_new(i,:) = getInvariants(C_mat(:,:,i),a0);
    
end
genPoints = genPoints_new;

%% Define inline function for rotation matrices
Rx = @(tx) [1, 0, 0;
    0, cos(tx), - sin(tx);
    0, sin(tx),   cos(tx)];

Ry = @(ty)[cos(ty),  0,   sin(ty);
    0,  1,   0;
    -sin(ty), 0,   cos(ty)];

Rz =@(tz) [cos(tz), - sin(tz), 0;
    sin(tz), cos(tz),   0;
    0,         0,       1];


%% Algorithm 3 of paper
steps = 20001;
maxStepSize =  2*pi;
a = 0;
b = 2*pi;
% Sample random initial angles
angles = (b-a).*rand(size(genPoints,1),3)  + a;
for i=1:steps
    % Display progress
    if mod(i,2000)==1
        ST = ['Step: ', num2str(i), ' of ', num2str(steps)];
        disp(ST);
    end
    
    % an array of d normally distributed random variables using Box-Muller
    % transform (eq. 50)
    u = normrnd(0,1,size(genPoints,1),3);
    r = sum(u.^2,2) .^(0.5);
    n = u./r;
    
    % Find nearest neighbors of current pseudo-invariants
    [~, dist] = knnsearch(genPoints(:,4:5),genPoints(:,4:5),'k',2);
    currentDistance = dist(:,2);
    
    % Sample stepsize
    stepsize = 0 + (maxStepSize-0).*rand(size(genPoints,1),1);
    
    % Next possible rotation angles
    possibleAngles = angles + stepsize.*n;
    
    % Get new possible pseudo invariant values
    C_test = zeros(3,3,size(genPoints,1));
    possiblePosition = zeros(size(genPoints,1),5);
    for jj=1:size(genPoints,1)
        C_test(:,:,jj)  = Rx(possibleAngles(jj,1))'*Ry(possibleAngles(jj,2))'*Rz(possibleAngles(jj,3))'  * C_mat(:,:,jj) *Rz(possibleAngles(jj,3))*Ry(possibleAngles(jj,2))*Rx(possibleAngles(jj,1))  ;
        possiblePosition(jj,:) = getInvariants(C_test(:,:,jj),a0);
    end
    
    % Check if distance to nearest neighbor increased
    [~, dist] = knnsearch(genPoints(:,4:5),possiblePosition(:,4:5),'k',1);
    possibleDistance = dist(:,1);
    test1 = possibleDistance>currentDistance;
    
    % Update angles and positions
    angles = angles + stepsize.*n.*test1;
    idx = find(test1>0);
    genPoints(idx,:) = possiblePosition(idx,:);
    
    % Update stepsize
    maxStepSize = 0.99995*maxStepSize;
end



end

