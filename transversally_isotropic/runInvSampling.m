%% 
% Code for sampling evenly spaced samples in transverally isotropic invariant space 
% using simulated annealing
% For more information refer to https://arxiv.org/pdf/2109.11028.pdf
% Fuhg, Jan Niklas, and Nikolaos Bouklas. 
% "On physics-informed data-driven isotropic and anisotropic constitutive 
% models through probabilistic machine learning and space-filling sampling." 
% arXiv preprint arXiv:2109.11028 (2021).
%
% Check Section 4 of the paper for more information

clear all;
rng(0);
%% As pointed out in the paper the process is modular, so we first create 
% evenly-spaced samples in invariant space using Algorithm 2 of the paper. 
% Then we use these samples to generate evenly-spaced combinations of 
% principal and pseudo-invariants needed for transversally isotropic materials.


% Define a preferred direction vector
a0 = [1;2;1];
a0 = a0/norm(a0);

%% First principal invariant space points
% Create the hull in which the sample points should lie
% Number of samples to generate hull 
% (The process can be significantly sped up by reducing this number)
n_points_hull = 25000;
% Delta of the deformation gradient as defined by eq.(45) in the paper
delta = 0.175;
% Get points that describe hull (compare Algorithm 1)
invHull = getConvHull(n_points_hull, delta,a0);


% Number of points to sample
n_points_sample = 50;

% Sampling algorithm for principal invariants (compare Algorithm 2)
disp('Sampling principal invariant space');
[genPoints] = FunSimAnnealing(n_points_sample,invHull);


%% Next using these evenly-spaced points generate space-filling points in 
% pseudo-invariant space using Algorithm 3
disp('Sampling pseudo invariant space');
[genPoints] = FunSimAnnealing45(genPoints,a0);


%% As an example you can obtain a right Cauchy-Green tensor from the set
% (I1,I2,I3,I4,I5) using the following function (eq. (55)). Here for the
% 5th point of the generated dataset (randomly chosen)
[C] = funGetCfromInv(genPoints(5,:),delta,a0);



%% Plotting and saving
figure
scatter3(invHull(:,1), -invHull(:,2), invHull(:,3),1, 'filled','b'); hold on;
scatter3(genPoints(2:n_points_sample,1), -genPoints(2:n_points_sample,2), genPoints(2:n_points_sample,3),35, 'filled', 'r'); hold on;
scatter3(genPoints(1,1), -genPoints(1,2), genPoints(1,3),50, 'filled', 'k'); hold on;
grid on;
xlabel('$I_{1}$', 'Interpreter', 'Latex','FontSize',18);
ylabel('$I_{2}$', 'Interpreter', 'Latex','FontSize',18);
zlabel('$I_{3}$', 'Interpreter', 'Latex','FontSize',18);
legend({'Convex hull', 'Sampled points', 'Undeformed configuration'}, 'Interpreter', 'Latex', 'Location', 'best')
set(gca,'FontSize',18)
%view(-125,24)
box on;
set(gcf,'Renderer','Painter' ,'Position',[50 50 850 600]);
STOpti1 = ['I1I3Init.pdf'];
saveas(gcf,STOpti1)
STCrop = ['pdfcrop ',STOpti1, ' ',STOpti1];
system(STCrop)

figure
scatter(invHull(:,4), invHull(:,5),1, 'filled','b'); hold on;
scatter(genPoints(2:n_points_sample,4), genPoints(2:n_points_sample,5),35, 'filled', 'r'); hold on;
scatter(genPoints(1,4), genPoints(1,5),50, 'filled', 'k'); hold on;
grid on;
xlabel('$I_{4}$', 'Interpreter', 'Latex','FontSize',18);
ylabel('$I_{5}$', 'Interpreter', 'Latex','FontSize',18);
legend({'Convex hull', 'Sampled points', 'Undeformed configuration'}, 'Interpreter', 'Latex', 'Location', 'best')
set(gca,'FontSize',18)
%view(-125,24)
box on;
set(gcf,'Renderer','Painter' ,'Position',[50 50 850 600]);
STOpti1 = ['I4I5Init.pdf'];
saveas(gcf,STOpti1)
STCrop = ['pdfcrop ',STOpti1, ' ',STOpti1];
system(STCrop)


ST = ['InvGenerated_', num2str(n_points_sample), '.mat'];
save(ST,'genPoints');