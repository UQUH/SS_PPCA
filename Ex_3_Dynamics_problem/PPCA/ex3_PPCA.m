%% Apply SS-PPCA method to the Linear Dynamics Problem
% The problem setup follows the reference paper.

seed = 42; % Set random seed for reproducibility
rng(seed);

%% Load System Matrices and External Force
addpath('~/OneDrive/laptop/PhD/SROB/SS_PPCA/Ex_3_Dynamics_problem/LS_Dyna_model')

load("system_matrices.mat") % Load Mass, Stiffness, and Damping matrices
dof = 42342; % Degrees of freedom
beta_coff = 6.366e-06; % Raleigh damping coffecient (proportional to stiffness)

load("force.mat") % Load external force input

%% Time Discretization
N = 1537; % Number of time steps
dt = 5e-5; % Time step size (s)
t = 1000 * dt * (0:N-1); % Time vector (converted to ms)

%% Load High-Dimensional Model (HDM) Data & Compute Reduced Order Basis
load("displacement.mat",'U','Udot') % Load displacement & velocity data

% Perform Singular Value Decomposition (SVD)
[V,D,~] = svd(U, "econ"); % Compute economy-sized SVD of displacement matrix

num_modes = 150; % Number of dominant modes retained
V_r = V(:,1:num_modes); % Reduced-order basis

% Project System Matrices onto the Subspace
M_r = V_r' * M * V_r; % Reduced Mass matrix
K_r = V_r' * K * V_r; % Reduced Stiffness matrix
F_r = V_r' * F; % Reduced Force vector
D_r = D(1:num_modes,1:num_modes); % Truncated singular values

% Extract Global Reduced Basis
k = 10; % Number of modes for global basis
V_global = V(:,1:k); % Reduced-order basis for global analysis

%% load Reduced-order Model (ROM) data
load("displacement_rom_k_10.mat","v_ROM")

%% Specify Quantities of Interest (QoI)
rand_dof = 530; % Randomly selected degree of freedom (DOF)
dof_obs = 23374; % Observed degree of freedom (DOF) for analysis

qoi_HDM = Udot(dof_obs,:); % Velocity at observed DOF from HDM
qoi_HDM_rand = Udot(rand_dof,:); % Velocity at random DOF from HDM

qoi_ROM = v_ROM(dof_obs,:); % Velocity at observed DOF from ROM
qoi_ROM_rand = v_ROM(rand_dof,:); % Velocity at random DOF from ROM

clear U Udot Uddot v_ROM % Clean Up Unused Variables

%% Optimization for Hyperparameter Selection
dist_HDM_ROM = sqrt(sum((qoi_HDM - qoi_ROM).^2, "all")); % Compute L2 distance between HDM & ROM QoI

tic; % Start timing

betaEval = []; % Initialize empty array for beta evaluations
boundBeta = [1, 100]; % Define search bounds for beta optimization

options = optimset('TolX', 0.1); % Set optimization options (tolerance for fminbnd)

% Optimize beta using fminbnd
[betaOpt, fval] = fminbnd(@(beta) fMC(beta, V_r, D_r, K_r, M_r, ...
    beta_coff, num_modes, F_r, qoi_ROM, dist_HDM_ROM, k), ...
    boundBeta(1), boundBeta(2), options);

toc; % End timing

beta_opt = floor(betaOpt); % Round down optimal beta value (beta_opt = 39)

%% Compute Stochastic Reduced-Order Basis (SROB) and Solve System
store_node_states_SROM = zeros(18, N, 1000); % Store node states (displacement, velocity, acceleration)
store_qoi_SROM = zeros(1000, N); % Store QoI for observed DOF (SROM)
store_qoi_rand_SROM = zeros(1000, N); % Store QoI for random DOF (SROM)

for i = 1:1000
    Z = randn(num_rank_U, beta_opt);  % Generate random perturbation matrix
    mat = D_r * Z;
    
    [W, ~, ~] = svds(mat, k); % Perform truncated SVD
    
    M_SROM = W' * M_r * W; % Reduced mass matrix
    K_SROM = W' * K_r * W; % Reduced stiffness matrix
    C_SROM = beta_coff * K_SROM; % Reduced damping matrix
    F_SROM = W' * F_r; % Reduced force vector
    
    Qo = zeros(k,1); % Initial conditions
    Qdoto = zeros(k,1);
    
    % Solve system using Newmark-beta time integration scheme
    [Q, Qdot, Qddot] = Newmark(M_SROM, C_SROM, K_SROM, F_SROM, dt, N, Qo, Qdoto);
    
    % Project solutions back to full space
    y_SROM = V_r * W * Q; % Displacement
    v_SROM = V_r * W * Qdot; % Velocity
    a_SROM = V_r * W * Qddot; % Acceleration
    
    % Store QoI for observed and random DOFs
    store_qoi_rand_SROM(i, :) = v_SROM(rand_dof, :);
    store_qoi_SROM(i, :) = v_SROM(dof_obs, :);
    
    % Store state variables for selected DOFs
    store_node_states_SROM(1:6, :, i) = y_SROM(dof_obs:dof_obs+5, :); % Displacement
    store_node_states_SROM(7:12, :, i) = v_SROM(dof_obs:dof_obs+5, :); % Velocity
    store_node_states_SROM(13:18, :, i) = a_SROM(dof_obs:dof_obs+5, :); % Acceleration
end

% Clean up temporary variables
clear y_SROM v_SROM a_SROM

%% Plot: Velocity at Observed DOF with 95% Prediction Interval
pc = 0.95; % Confidence level
alpha = 1 - pc;

vel_lower_PPCA = quantile(store_qoi_SROM, alpha / 2);% Compute 95% prediction interval
vel_upper_PPCA = quantile(store_qoi_SROM, 1 - alpha / 2);

mean_qoi_SROM = mean(store_qoi_SROM); %mean of SROM qoi
indicate_dof_type = 1; %indicator for velocity dof to help in label writing

% Generate the plot
f = custom_plot(vel_lower_PPCA, vel_upper_PPCA, qoi_HDM, qoi_ROM, mean_qoi_SROM, t,indicate_dof_type);

% Save figure in multiple formats
filename = "95_percent_prediction_PPCA_vel";
saveas(f, filename, 'fig'); % Save as MATLAB figure
saveas(f, filename, 'svg'); % Save as SVG format
exportgraphics(f, '95_percent_prediction_PPCA_vel.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Calculate the coverage and width for velocity
coverage_vel = coverage(vel_upper_PPCA,vel_lower_PPCA,qoi_HDM);%(after first initial 200 timesteps)
width_vel = mean(vel_upper_PPCA - vel_lower_PPCA);

%% Plot: Acceleration at Observed DOF with 95% Prediction Interval
acc_dof = 13; % DOF corresponding to acceleration

% Extract acceleration data and normalize
qoi_SROM_acc = store_node_states_SROM(acc_dof, :, :) / 1e4;
qoi_ROM_acc = store_node_states_ROM(acc_dof, :) / 1e4;
qoi_HDM_acc = store_node_states_HDM(acc_dof, :) / 1e4;

mean_SROM_acc = mean(qoi_SROM_acc, 3); % Compute mean acceleration from SROM

pc = 0.95;% Set confidence level
alpha = 1 - pc;

reshape_qoi_SROM_acc = squeeze(qoi_SROM_acc)';% Reshape SROM data for quantile computation

% Compute 95% prediction interval
acc_lower_PPCA = quantile(reshape_qoi_SROM_acc, alpha / 2);
acc_upper_PPCA = quantile(reshape_qoi_SROM_acc, 1 - alpha / 2);

% Indicator for acceleration DOF to assist in label writing
indicate_dof_type = 2;

% Generate plot
f = custom_plot(acc_lower_PPCA, acc_upper_PPCA, qoi_HDM_acc, ...
    qoi_ROM_acc, mean_SROM_acc, t, indicate_dof_type);

% Save figure in multiple formats
filename = "95_percent_prediction_PPCA_acc";
saveas(f, filename, 'fig'); % Save as MATLAB figure
saveas(f, filename, 'svg'); % Save as SVG format
exportgraphics(f, '95_percent_prediction_PPCA_acc.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Calculate the coverage and width for acceleration 
coverage_acc = coverage(acc_upper_PPCA,acc_lower_PPCA,qoi_HDM_acc); %(after first initial 200 timesteps)
width_acc = mean(acc_upper_PPCA - acc_lower_PPCA);

%% Plot: Displacement at Observed DOF with 95% Prediction Interval
disp_dof = 1; % DOF corresponding to displacement

% Extract displacement data and scale
qoi_SROM_disp = store_node_states_SROM(disp_dof, :, :) * 1e3;
qoi_ROM_disp = store_node_states_ROM(disp_dof, :) * 1e3;
qoi_HDM_disp = store_node_states_HDM(disp_dof, :) * 1e3;

mean_qoi_SROM_disp = mean(qoi_SROM_disp, 3);% Compute mean displacement from SROM

pc = 0.95;% Set confidence level
alpha = 1 - pc;

reshape_qoi_SROM_disp = squeeze(qoi_SROM_disp)';% Reshape SROM data for quantile computation

% Compute 95% prediction interval
disp_lower_PPCA = quantile(reshape_qoi_SROM_disp, alpha / 2);
disp_upper_PPCA = quantile(reshape_qoi_SROM_disp, 1 - alpha / 2);

% Indicator for displacement DOF to assist in label writing
indicate_dof_type = 3;

% Generate the plot
f = custom_plot(disp_lower_PPCA, disp_upper_PPCA, qoi_HDM_disp, ...
    qoi_ROM_disp, mean_qoi_SROM_disp, t, indicate_dof_type);

% Save figure in multiple formats
filename = "95_percent_prediction_PPCA_disp";
saveas(f, filename, 'svg'); % Save as SVG format
saveas(f, filename, 'fig'); % Save as MATLAB figure
exportgraphics(f, '95_percent_prediction_PPCA_disp.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Calculate the coverage and width for displacment
coverage_disp = coverage(disp_upper_PPCA,disp_lower_PPCA,qoi_HDM_disp);%(after first initial 200 timesteps)
width_disp = mean(disp_upper_PPCA - disp_lower_PPCA);

%% Plot: Velocity at a Random DOF with 95% Prediction Interval
mean_qoi_SROM_rand = mean(store_qoi_rand_SROM, 1); % Compute mean velocity for random DOF

% Compute 95% prediction interval
y_lower_rand = quantile(store_qoi_rand_SROM, alpha / 2);
y_upper_rand = quantile(store_qoi_rand_SROM, 1 - alpha / 2);

% Generate the plot
f = custom_plot(y_lower_rand, y_upper_rand, qoi_HDM_rand, ...
    qoi_ROM_rand, mean_qoi_SROM_rand, t);

% Save figure in multiple formats
filename = "95_percent_prediction_PPCA_rand_dof";
saveas(f, filename, 'fig'); % Save as MATLAB figure
saveas(f, filename, 'svg'); % Save as SVG format
exportgraphics(f, '95_percent_prediction_PPCA_rand_dof.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Compare Width of NPM and PPCA

% Load NPM data for caculating width
addpath("~/OneDrive/laptop/PhD/SROB/SS_PPCA/Ex_3_Dynamics_problem/NPM")
load("width_data_NPM.mat"); % Load NPM data

% Compute difference between upper and lower bound for NPM
vel_diff_NPM = vel_upper_NPM - vel_lower_NPM;% Velocity
acc_diff_NPM = acc_upper_NPM - acc_lower_NPM; % Acceleration
disp_diff_NPM = disp_upper_NPM - disp_lower_NPM; %Displacement

% Compute difference between upper and lower bound for PPCA
vel_diff_PPCA = vel_upper_PPCA - vel_lower_PPCA;
acc_diff_PPCA = acc_upper_PPCA - acc_lower_PPCA;
disp_diff_PPCA = disp_upper_PPCA - disp_lower_PPCA;

% Compute Width Ratios (NPM/PPCA)
width_ratio_vel = mean(vel_diff_NPM ./ vel_diff_PPCA); % Velocity width ratio
width_ratio_acc = mean(acc_diff_NPM ./ acc_diff_PPCA); % Acceleration width ratio
width_ratio_disp = mean(disp_diff_NPM ./ disp_diff_PPCA); % Displacement width ratio

% Display Results
disp(['Width Ratio (Velocity): ', num2str(width_ratio_vel)]);
disp(['Width Ratio (Acceleration): ', num2str(width_ratio_acc)]);
disp(['Width Ratio (Displacement): ', num2str(width_ratio_disp)]);