function [genPoints] = FunSimAnnealing(n_points,hullPoints)
% inhull: Use of simulated annealing to generate space-filling points in
%         isotropic invariant space (Algorithm 2 of paper)
% 
% arguments: (input)
%  n_points: Number of points we want to generate with the simulated
%            annealing algorithm
%  hullPoints: Prior sampled points that define the convex hull of the
%              deformation gradient domain we aim to sample in (Algorithm 1)
%
% arguments: (output)
%  genPoints: Sampled points through the algorithm


% From the previously sampled points that define the convex hull 
% we look for n_points number of points that adhere to the physical
% restrictions (not strictly necessary anymore but a good initial check)
notFound = 1;
while notFound
    % Pick n_points random points of hullPoints
    r = randi([1 size(hullPoints,1)],n_points,1);
    tespPoints = hullPoints(r,:);
    
    % Check if physics constraints satisfied (eq (46)-(48))
    H = (1/9)*(tespPoints(:,1).^2 - 3*tespPoints(:,2));
    G = (1/3)*tespPoints(:,1).*tespPoints(:,2) - tespPoints(:,3) - 2*(tespPoints(:,1)./3).^3;
    beta = acos(-G./(2*(H).^(3/2)));
    lam31_sq = tespPoints(:,1)./3 + 2*(H).^(1/2).*cos(beta./3);
    lam21_sq = tespPoints(:,1)./3 - 2*(H).^(1/2).*cos((pi+beta)./3);
    lam11_sq = tespPoints(:,1)./3 - 2*(H).^(1/2).*cos((pi-beta)./3);
    test_physical = (imag(lam31_sq)==0) & (imag(lam21_sq)==0) & (imag(lam11_sq)==0);
    
    % If satisfied exit while loop
    if numel(find(test_physical==0))==0
        notFound = 0;
    end
end


% We want to enforce that one point of the newly generated points
% corredponds to the invariants of the undeformed configuration.
% So we replace the first point of the current set of point with the
% undeformed configuration point
undeformed = [3,3,1];
genPoints = tespPoints;
genPoints(1,:) = undeformed;


% Algorithm 2 of paper
steps = 10001;
maxStepSize =  1.0;
for i=1:steps
    % Display progress
    if mod(i,2000)==1
        ST = ['Step: ', num2str(i), ' of ', num2str(steps)];
        disp(ST);
    end
    
    % Generate random step directions
    rvals = 2*rand(size(genPoints,1),1)-1;
    elevation = asin(rvals);
    azimuth = 2*pi*rand(size(genPoints,1),1);
    [x,y,z] = sph2cart(azimuth,elevation,1.0);
    n = [x,y,z];
    
    % Find current nearest distances for each point of the set
    [~, dist] = knnsearch(genPoints,genPoints,'k',2);
    currentDistance = dist(:,2);
    
    % Sample step size
    stepsize = 0 + (maxStepSize-0).*rand(size(genPoints,1),1);
    
    % To leave undeformed configuration unchanged set stepsize(1)=0;
    stepsize(1)=0;
    
    % Possible next position of points
    possiblePosition = genPoints + stepsize.*n;
    
    % Check if new positions are inside convex hull
    in = inhull(possiblePosition,hullPoints);
    
    % Check if new positions increase the distance to nearest neighbore
    [~, dist] = knnsearch(genPoints,possiblePosition,'k',2);
    possibleDistance = dist(:,2);
    test1 = possibleDistance>currentDistance;
    
    % Check for physicality of new positions
    H = (1/9)*(possiblePosition(:,1).^2 - 3*possiblePosition(:,2));
    G = (1/3)*possiblePosition(:,1).*possiblePosition(:,2) - possiblePosition(:,3) - 2*(possiblePosition(:,1)./3).^3;
    beta = acos(-G./(2*(H).^(3/2)));
    lam31_sq = possiblePosition(:,1)./3 + 2*(H).^(1/2).*cos(beta./3);
    lam21_sq = possiblePosition(:,1)./3 - 2*(H).^(1/2).*cos((pi+beta)./3);
    lam11_sq = possiblePosition(:,1)./3 - 2*(H).^(1/2).*cos((pi-beta)./3);
    test2 = (imag(lam31_sq)==0) & (imag(lam21_sq)==0) & (imag(lam11_sq)==0);
    
    % New point positions based on checks from above
    genPoints = genPoints + stepsize.*n.*in.*test1.*test2;
    
    % Reduce the step size
    maxStepSize = 0.9995*maxStepSize;
    
end


end

