%% Parametric Linear static problem: SS-PPCA
rng(42);% Set random seed for reproducibility

%% Add Model folder into path
% Get the parent directory of the current working directory
parentDir = fileparts(pwd);

% Define the target folder (which exists in the parent directory)
targetFolder = fullfile(parentDir, 'Model');

% Check if the directory exists before adding it to the path
if isfolder(targetFolder)
    addpath(targetFolder);
    disp(['Successfully added path: ', targetFolder]);
else
    warning(['The directory "', targetFolder, '" does not exist.']);
end

%% Generate Data using the Model
[K,K_bc,V_global,f_vec,phi_part,y_HDM,y_ROM,n,k,...
    zeta,num_models] = Ex1Model;

%% Perfrom SVD on perturbed HDM Models data
[V_star,D,~] = svd(y_HDM,"econ"); %economy-sized SVD
D = diag(D); % transform into column vector

%% Optimization for Hyperparameter
% Compute L2 distance between HDM & ROM Observed values
dist_HDM_ROM = sqrt(sum((y_HDM - y_ROM).^2,"all")); 

tic; % Start timing

betaEval = []; % Initialize empty array for beta evaluations
boundBeta = [k, 50]; % Define search bounds for beta optimization

options = optimset('TolX', 0.1); % Set optimization options (tolerance for fminbnd)

% Optimize beta using fminbnd
[betaOpt, fval] = fminbnd(@(beta) fMC(beta,V_star,D,K,f_vec,...
    dist_HDM_ROM,y_ROM,num_models,n,k), ...
    boundBeta(1), boundBeta(2), options);

% Compute and display execution time
Opt_time_PPCA = toc; % End timing
fprintf('Hyper paramaeter optimization completed in %.2f seconds.\n', Opt_time_PPCA);
beta_opt = floor(betaOpt); % Round down optimal beta value (beta_opt = )

%% Performance of the mu-parametric SROM at test parameter value
mu_test = [0.5,0.5,0.5,0.5,1];
g_test = sum(mu_test .* phi_part,2);
g_max_test = max(abs(g_test));
f_vec_test = (1/g_max_test)*g_test;

% Solve HDM
y_HDM_test = zeros(n,1);
y_HDM_test(2:end-1) = K_bc\f_vec_test(2:end-1); %HDM

% Solve ROM
K_ROM = V_global' * K * V_global;             % Reduced stiffness matrix
q_test = K_ROM \ (V_global' * f_vec_test);     % Solve for reduced coefficients
y_ROM_test = V_global * q_test;                   % Map back to full space

%% Plot properties
set(0, 'DefaultAxesFontSize', 12);                    % Set default axes font size
set(0, 'DefaultAxesFontName', 'Times');               % Set default font to Times
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');   % Use LaTeX for axis tick labels
set(0, 'DefaultLegendInterpreter', 'latex');          % Use LaTeX for legend text
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % Redundant but ensures LaTeX ticks

width_plot = 540;  
height_plot = width_plot / 3; 
th = 1.4; 

%% To store figures
% Create the folder if it doesn't exist
if ~exist(fullfile(pwd, 'Results'), 'dir')
    mkdir(fullfile(pwd, 'Results'));
end

%% Plot: HDM vs ROM at test parameters
figure('Units', 'points', 'Position', [1, 1, width_plot, height_plot]);  % Set figure size
plot(zeta, y_HDM_test, 'k', 'LineWidth', th);                % Plot HDM in black
hold on;
p2 = plot(zeta, y_ROM_test, 'LineWidth', th);                % Plot ROM 
p2.Color = '#ee3a2b';                                      % Set ROM line color using hex code

xlabel('X');                
ylabel('Displacement', Interpreter='latex');
ylim([-1.2 2.2]*1e-3);   

legend('HDM', 'ROM', 'Box', 'off', 'Location', 'northeast'); % Add legend without box
filename = "HDM_ROM_Ex2.pdf";  
output_path = fullfile(pwd, 'Results', filename);% Combine into a full relative path
exportgraphics(gcf, output_path, 'ContentType', 'vector');% Export the figure

%% Generate stochastic reduced-order basis and solve SROM
v_sim = 1000;                               % Number of stochastic simulations (Monte Carlo samples)
y_SROM_test = zeros(n, v_sim);               % Preallocate matrix to store full-field solutions

for i = 1:v_sim
    Z = randn(num_models, beta_opt);           % Generate standard normal random matrix (num_models × beta_opt)
    mat = D.* Z;                               % Apply transformation
    [W, ~, ~] = svds(mat, k);                  % Compute top-k left singular vectors (low-rank basis)
    
    W = V_star * W;                            % Map reduced basis back to full space using V_star

    K_SROM = W' * K * W;                       % Project stiffness matrix to reduced-order space
    q_test = K_SROM \ (W' * f_vec_test);     % Solve for reduced coordinates (ROM solution)
    y_SROM_test(:, i) = W * q_test;             % Lift reduced solution back to full space
end

mean_y_srom_test = mean(y_SROM_test,2);

%% Plot: SROM Prediction with 95% Prediction Interval  
set(0,'DefaultAxesFontSize',12)
width_plot = 540;
height_plot = width_plot/3;

th = 1.4; % thickness of the line

pc = 0.95;
alpha = 1 - pc;
y_lower_PPCA = quantile(y_SROM_test',alpha/2);
y_upper_PPCA = quantile(y_SROM_test',1-alpha/2);

figure('units','points','position',[0,0,width_plot,height_plot])
h = fill([zeta', flip(zeta')], [y_lower_PPCA, flip(y_upper_PPCA)],'c','FaceAlpha',1);  % plot filled area
h.FaceColor = '#a6cce3';
h.EdgeColor = "none";
hold on 

p1 = plot(zeta,y_HDM_test,'k','LineWidth',th);

p2 = plot(zeta,y_ROM_test,'LineWidth',th);
p2.Color = '#ee3a2b';
p3 = plot(zeta,mean_y_srom_test,'LineWidth',th);
p3.Color = "#1e78b3"; 

xlabel('X')
xticks(0:0.2:1)
ylim([-1.2 2.2]*1e-3)
ylabel('Displacement');
legend([p1,p2,p3,h],{'HDM','ROM','SROM mean','SROM $95\%$ PI'},Location='northeast',Box='off',Interpreter='latex')

filename = "Prediction_PPCA_Ex2.pdf";  
output_path = fullfile(pwd, 'Results', filename);% Combine into a full relative path
exportgraphics(gcf, output_path, 'ContentType', 'vector');% Export the figure

%% Calculate the coverage and width for displacment
coverage_PPCA = coverage(y_upper_PPCA,y_lower_PPCA,y_HDM_test');
width_PPCA = mean(y_upper_PPCA - y_lower_PPCA);
