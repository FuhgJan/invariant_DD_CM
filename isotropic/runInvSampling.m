%% 
% Code for sampling evenly spaced samples in isotropic invariant space 
% using simulated annealing
% For more information refer to https://arxiv.org/pdf/2109.11028.pdf
% Fuhg, Jan Niklas, and Nikolaos Bouklas. 
% "On physics-informed data-driven isotropic and anisotropic constitutive 
% models through probabilistic machine learning and space-filling sampling." 
% arXiv preprint arXiv:2109.11028 (2021).
%
% Check Section 4 of the paper for more information

clear all;


%% Create the hull in which the sample points should lie
% Number of samples to generate hull
n_points_hull = 50000;
% Delta of the deformation gradient as defined by eq.(45) in the paper
delta = 0.175;
% Get points that describe hull (compare Algorithm 1)
invHull = getConvHull(n_points_hull, delta);


%% Using the hull generate space-filling sampling points inside it that
% that are physically possible

% Number of points to sample
n_points_sample = 200;

% Sampling algorithm for isotropic material symmetry (compare Algorithm 2)
[genPoints] = FunSimAnnealing(n_points_sample,invHull);



%% As an example you can obtain a right Cauchy-Green tensor from the triple
% (I1,I2,I3) using the following function (eq. (48), (49)). Here for the
% 5th point of the generated dataset (randomly chosen)
disp('Sampling principal invariant space')
[C] = getCFromInv(genPoints(5,1),genPoints(5,2),genPoints(5,3));


%% Plotting and saving
figure
scatter3(invHull(:,1), -invHull(:,2), invHull(:,3),1, 'filled','b'); hold on;
scatter3(genPoints(:,1), -genPoints(:,2), genPoints(:,3),35, 'filled', 'r'); hold on;
scatter3(genPoints(1,1), -genPoints(1,2), genPoints(1,3),50, 'filled', 'k'); hold on;
grid on;
xlabel('$I_{1}$', 'Interpreter', 'Latex','FontSize',18);
ylabel('$I_{2}$', 'Interpreter', 'Latex','FontSize',18);
zlabel('$I_{3}$', 'Interpreter', 'Latex','FontSize',18);
legend({'Convex hull', 'Sampled points'}, 'Interpreter', 'Latex', 'Location', 'best')
set(gca,'FontSize',18)
%view(-125,24)
box on;
set(gcf,'Renderer','Painter' ,'Position',[50 50 850 600]);


ST = ['InvGenerated_', num2str(n_points_sample), '.mat'];
save(ST,'genPoints');