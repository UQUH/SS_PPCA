%% Apply NPM method to the Linear Dynamics Problem
%%
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

%%
Node_data = readmatrix('Node_Data_0001_001');
Node_data(1,:) = []; %remove header
Node_data_spatial = Node_data(9998:end,1:4);
Node_data_spatial(9932:9991,:) = [];
Node_data_dof = Node_data(1:9997,2:end); %node data with dof
Node_data_dof(9932:9991,:) = [];

nodes = [3085, 3103, 3106, 3109, 3112, 3115, 3118, 3121,...
    3124, 3127, 3130, 3133, 3136, 3139, 3142, 3145,...
    3067, 3073, 3076, 3079, 3082, 3088, 3091, 3094,...
    3097, 3100, 3025, 3028, 3031, 3034, 3037, 3040,...
    3043, 3046, 3049, 3052, 3055, 3058, 3061, 3064,...
    3070, 3004, 3007, 3010, 3013, 3016, 3019, 3022];
for i = 1:48
    k = nodes(1,i);
    Node_data_dof(k,4:6) = 0;
end
N_0 = 9937;
Node_data_dof_vec = reshape(Node_data_dof,[N_0*6,1]);
[sorted_data,sorted_index] = sort(Node_data_dof_vec);
sorted_index = sorted_index(sorted_data~=0);
sorted_data(sorted_data == 0) = [];

%% Load High-Dimensional Model (HDM) Data & Compute Reduced Order Basis
load("displacement.mat",'U','Udot') % Load displacement & velocity data

% Perform Singular Value Decomposition (SVD)
[V,D,~] = svd(U, "econ"); % Compute economy-sized SVD of displacement matrix

% Extract Global Reduced Basis
k = 10; % Number of modes for global basis
V_global = V(:,1:k); % Reduced-order basis for global analysis

%% load Reduced-order Model (ROM) data
load("displacement_rom_k_10.mat","v_ROM")

%% Specify Quantities of Interest (QoI)
dof_obs = 23374; % Observed degree of freedom (DOF) for analysis

qoi_HDM = Udot(dof_obs,:); % Velocity at observed DOF from HDM
qoi_ROM = v_ROM(dof_obs,:); % Velocity at observed DOF from ROM

clear U Udot Uddot v_ROM % Clean Up Unused Variables

%% Pre computation as per the reference paper Soize and Farhat (2017)
d = 3;
L = zeros(d,1);
for i = 1:d
    L(i) = max(Node_data_spatial(:,i+1)) - min(Node_data_spatial(:,i+1));
end

vp = 20;
k_v = zeros(1,vp);
lambda_v = zeros(1,vp);
for i = 1:vp
    k_v(i) = -1 + (i - 0.5)*(2/vp);
    if k_v(i) == 0
        lambda_v(i) = (2/vp)*1;
    else
        lambda_v(i) = (2/vp)*(1-abs(k_v(i)));
    end
end
N_0 = 9937;
x = Node_data_spatial(:,2:end);

%% The optimization problem have been solved and the optimal
% hyper-parameter have been saved to reduce the computational time while
% reproducing the results. If someone wants to retrain then uncomment the
% code below and retrain. Please note that it might take a day or more to
% retrain !!
%% Solving the Optimization Problem - Stage 1
% tic; % Start timing
% 
% % Initial parameter values
% s_0 = 0.017; % Initial s value
% beta_0 = 0.2; % Initial beta value
% sigma_0 = ones(1, k); % Initial sigma values
% alpha_0 = [s_0, beta_0, sigma_0]; % Initial parameter vector
% 
% % Define objective function for Stage 1
% fun = @(alpha) objective_fun_NPM_stage1(alpha, M, K, beta_coff, ...
%     sorted_index, F, qoi_HDM, qoi_ROM, V_global, lambda_v, k_v, L, x);
% 
% % Define lower and upper bounds for optimization
% lb = [0.01, 0.15, 0.01 * ones(1, k)]; % Lower bounds
% ub = [1, 0.3, 20 * ones(1, k)]; % Upper bounds
% 
% % Set optimization options
% opts = optimoptions('fmincon', 'Display', 'iter');
% 
% % Create optimization problem
% problem = createOptimProblem('fmincon', 'objective', fun, ...
%     'x0', alpha_0, 'lb', lb, 'ub', ub, 'options', opts);
% 
% % Solve optimization problem
% alpha_stage_1 = fmincon(problem);
% 
% % Save optimized parameters
% save("alpha_opt_stage1.mat", "alpha_stage_1");
% 
% %% Stage 2 - Solving the Optimization Problem
% % Initialize optimization variables
% alpha_0 = zeros(1, k * (k - 1) * 0.5); % Initial parameter vector
% 
% % Define objective function for Stage 2
% fun = @(alpha) objective_fun_NPM_stage2(alpha, M, K, beta_coff, ...
%     sorted_index, F, qoi_HDM, qoi_ROM, V_global, lambda_v, k_v, L, x);
% 
% % Define lower and upper bounds (empty if no constraints)
% lb = [];
% ub = [];
% 
% % Set optimization options
% opts = optimoptions('fmincon', 'Display', 'iter');
% 
% % Create optimization problem
% problem = createOptimProblem('fmincon', 'objective', fun, ...
%     'x0', alpha_0, 'lb', lb, 'ub', ub, 'options', opts);
% 
% % Solve optimization problem
% alpha_stage_2 = fmincon(problem);
% 
% % Save optimized parameters
% save("alpha_opt_stage2.mat", "alpha_stage_2");
% 
% %% Stage 3 - Solving the Optimization Problem
% % Initialize optimization variable
% beta_0 = alpha_stage_1(2); % Extract beta from Stage 1 result
% alpha_0 = beta_0; % Initial guess for optimization
% 
% % Define objective function for Stage 3
% fun = @(alpha) objective_fun_NPM_stage3(alpha, M, K, beta_coff, ...
%     sorted_index, F, qoi_HDM, qoi_ROM, V_global, lambda_v, k_v, L, x);
% 
% % Define lower and upper bounds
% lb = 0.15;
% ub = 0.3;
% 
% % Set optimization options
% opts = optimoptions('fmincon', 'Display', 'iter');
% 
% % Create optimization problem
% problem = createOptimProblem('fmincon', 'objective', fun, ...
%     'x0', alpha_0, 'lb', lb, 'ub', ub, 'options', opts);
% 
% % Solve optimization problem
% alpha_stage_3 = fmincon(problem);
% 
% % Save optimized parameters
% save("alpha_opt_stage3.mat", "alpha_stage_3");
% 
% %% Stage 4 - Solving the Optimization Problem
% tic; % Start timing
% 
% % Initialize optimization variables
% s_0 = alpha_stage_1(1); % Extract s from Stage 1 result
% sigma_0 = [alpha_stage_1(3:end), alpha_stage_2]; % Combine sigma values from previous stages
% alpha_0 = [s_0, sigma_0]; % Initial parameter vector
% 
% % Define objective function for Stage 4
% fun = @(alpha) objective_fun_NPM_stage4(alpha, M, K, beta_coff, ...
%     sorted_index, F, qoi_HDM, qoi_ROM, V_global, lambda_v, k_v, L, x);
% 
% % Define lower and upper bounds
% lb = [0.01, 0.01 * ones(1, k)];
% ub = [1, 20 * ones(1, k)];
% 
% % Set optimization options
% opts = optimoptions('fmincon', 'Display', 'iter');
% 
% % Create optimization problem
% problem = createOptimProblem('fmincon', 'objective', fun, ...
%     'x0', alpha_0, 'lb', lb, 'ub', ub, 'options', opts);
% 
% % Solve optimization problem
% alpha_stage_4 = fmincon(problem);
% 
% % Save optimized parameters
% save("alpha_opt_stage4.mat", "alpha_stage_4");
% 
% toc; % End timing

%% Load Optimized Parameters from All Stages
load("alpha_opt_stage1.mat", "alpha_stage_1"); % Stage 1 results
load("alpha_opt_stage2.mat", "alpha_stage_2"); % Stage 2 results
load("alpha_opt_stage3.mat", "alpha_stage_3"); % Stage 3 results
load("alpha_opt_stage4.mat", "alpha_stage_4"); % Stage 4 results

% Construct Final Optimized Parameter Vector
alpha_opt = [alpha_stage_4(1), alpha_stage_3, alpha_stage_4(2:end)];

% Create an empty upper triangular matrix
alpha_opt_mat = triu(zeros(k)); 

% Populate the upper triangular matrix (including diagonal) using reshaping
alpha_opt_mat(triu(true(k))) = alpha_opt(3:end);

%% Performance of the Stochastic Reduced-Order Model (SROM)
v_sim = 1000; % Number of simulations
G = zeros(dof, k); % Initialize G matrix
G_node = zeros(N_0, 6); % Initialize node-level G matrix
m_dof = 6; % Number of degrees of freedom per node

% Storage for QoI and node states
store_qoi_SROM = zeros(v_sim, N);
store_node_states_SROM = zeros(18, N, v_sim);

for i = 1:v_sim
    for j = 1:k
        for l = 1:m_dof
            
            a_j_x = pi * k_v / (L(1) * alpha_opt(2));
            a_j_y = pi * k_v / (L(2) * alpha_opt(2));
            a_j_z = pi * k_v / (L(3) * alpha_opt(2));

            a_j_x_i = x(:, 1) .* a_j_x;
            a_j_y_i = x(:, 2) .* a_j_y;
            a_j_z_i = x(:, 3) .* a_j_z;

            G_x = cos(a_j_x_i);
            G_y = cos(a_j_y_i);
            G_z = cos(a_j_z_i);
            S_x = sin(a_j_x_i);
            S_y = sin(a_j_y_i);
            S_z = sin(a_j_z_i);

            C_S_vec = randn(vp, 2);
            C_S_mod = sqrt(C_S_vec(:, 1).^2 + C_S_vec(:, 2).^2);
            C_vec = C_S_vec(:, 1) ./ C_S_mod;
            S_vec = C_S_vec(:, 2) ./ C_S_mod;

            W_bar = sqrt(lambda_v)' .* abs(randn(vp, 1));
            W_C = W_bar .* C_vec;
            W_S = W_bar .* S_vec;

            g_X_i = G_x * W_C + S_x * W_S;
            g_Y_i = G_y * W_C + S_y * W_S;
            g_Z_i = G_z * W_C + S_z * W_S;

            G_node(:, l) = (1 / sqrt(dof)) * g_X_i .* g_Y_i .* g_Z_i;
        end
        
        G_reshaped = reshape(G_node, [N_0 * 6, 1]);
        G(:, j) = G_reshaped(sorted_index);
    end

    U_mat = G * alpha_opt_mat;
    %A = U_mat; A is same as U as there is no boundary conditions repalced
    
    D = (V_global' * U_mat + U_mat' * V_global) / 2;
    
    Z = U_mat - V_global * D;
    
    H_s_Z = (eye(k) + (alpha_opt(1)^2) * (Z' * Z))^(-0.5);
    
    % Compute final reduced basis W
    W = (V_global + alpha_opt(1) * Z) * H_s_Z;

    % Compute reduced-order system matrices
    M_SROM = W' * M * W; % Reduced mass matrix
    K_SROM = W' * K * W; % Reduced stiffness matrix
    C_SROM = beta_coff * K_SROM; % Reduced damping matrix
    F_SROM = W' * F; % Reduced force vector

    % Initialize initial conditions
    Qo = zeros(k, 1); % Initial displacement
    Qdoto = zeros(k, 1); % Initial velocity

    % Solve the reduced system using the Newmark-beta method
    [Q, Qdot, Qddot] = Newmark(M_SROM, C_SROM, K_SROM, F_SROM, dt, N, Qo, Qdoto);

    % Project reduced solutions back to full space
    y_SROM = W * Q; % Displacement
    v_SROM = W * Qdot; % Velocity
    a_SROM = W * Qddot; % Acceleration

    % Store quantities of interest (QoI)
    store_qoi_SROM(i, :) = v_SROM(dof_obs, :);

    % Store node states (displacement, velocity, acceleration)
    store_node_states_SROM(1:6, :, i) = y_SROM(dof_obs:dof_obs+5, :); % Displacement
    store_node_states_SROM(7:12, :, i) = v_SROM(dof_obs:dof_obs+5, :); % Velocity
    store_node_states_SROM(13:18, :, i) = a_SROM(dof_obs:dof_obs+5, :); % Acceleration
end
clear y_SROM v_SROM a_SROM % Clean Up Unused Variables
 %% Plot: Velocity at Observed DOF with 95% Prediction Interval
pc = 0.95; % Confidence level
alpha = 1 - pc;

vel_lower_NPM = quantile(store_qoi_SROM, alpha / 2);% Compute 95% prediction interval
vel_upper_NPM = quantile(store_qoi_SROM, 1 - alpha / 2);

mean_qoi_SROM = mean(store_qoi_SROM); %mean of SROM qoi
indicate_dof_type = 1; %indicator for velocity dof to help in label writing

% Generate the plot
f = custom_plot(vel_lower_NPM, vel_upper_NPM, qoi_HDM, qoi_ROM, mean_qoi_SROM, t,indicate_dof_type);

% Save figure in multiple formats
filename = "95_percent_prediction_NPM_vel";
saveas(f, filename, 'fig'); % Save as MATLAB figure
saveas(f, filename, 'svg'); % Save as SVG format
exportgraphics(f, '95_percent_prediction_NPM_vel.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Calculate the coverage and width for velocity 
coverage_vel = coverage(vel_upper_NPM,vel_lower_NPM,qoi_HDM);%(after first initial 200 timesteps)
width_vel = mean(vel_upper_NPM - vel_lower_NPM);

%% Plot: Acceleration at Observed DOF with 95% Prediction Interval
acc_dof = 13; % DOF corresponding to acceleration

% Extract acceleration data and normalize
qoi_SROM_acc = store_node_states_SROM(acc_dof, :, :) / 1e4;
qoi_ROM_acc = store_node_states_ROM(acc_dof, :) / 1e4;
qoi_HDM_acc = store_node_states_HDM(acc_dof, :) / 1e4;

mean_SROM_acc = mean(qoi_SROM_acc, 3); % Compute mean acceleration from SROM

pc = 0.95;% Set confidence level
alpha = 1 - pc;

reshape_qoi_SROM_disp = squeeze(qoi_SROM_acc)';% Reshape SROM data for quantile computation

% Compute 95% prediction interval
acc_lower_NPM = quantile(reshape_qoi_SROM_disp, alpha / 2);
acc_upper_NPM = quantile(reshape_qoi_SROM_disp, 1 - alpha / 2);

% Indicator for acceleration DOF to assist in label writing
indicate_dof_type = 2;

% Generate plot
f = custom_plot(acc_lower_NPM, acc_upper_NPM, qoi_HDM_acc, ...
    qoi_ROM_acc, mean_SROM_acc, t, indicate_dof_type);

% Save figure in multiple formats
filename = "95_percent_prediction_NPM_acc";
saveas(f, filename, 'fig'); % Save as MATLAB figure
saveas(f, filename, 'svg'); % Save as SVG format
exportgraphics(f, '95_percent_prediction_NPM_acc.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Calculate the coverage and width for acceleration 
coverage_acc = coverage(acc_upper_NPM,acc_lower_NPM,qoi_HDM_acc);%(after first initial 200 timesteps)
width_acc = mean(acc_upper_NPM - acc_lower_NPM);

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
disp_lower_NPM = quantile(reshape_qoi_SROM_disp, alpha / 2);
disp_upper_NPM = quantile(reshape_qoi_SROM_disp, 1 - alpha / 2);

% Indicator for displacement DOF to assist in label writing
indicate_dof_type = 3;

% Generate the plot
f = custom_plot(disp_lower_NPM, disp_upper_NPM, qoi_HDM_disp, ...
    qoi_ROM_disp, mean_qoi_SROM_disp, t, indicate_dof_type);

% Save figure in multiple formats
filename = "95_percent_prediction_NPM_disp";
saveas(f, filename, 'svg'); % Save as SVG format
saveas(f, filename, 'fig'); % Save as MATLAB figure
exportgraphics(f, '95_percent_prediction_NPM_disp.pdf', 'ContentType', 'vector'); % Save as vector PDF

%% Calculate the coverage and width for displacment 
coverage_disp = coverage(disp_upper_NPM,disp_lower_NPM,qoi_HDM_disp);%(after first initial 200 timesteps)
width_disp = mean(disp_upper_NPM - disp_lower_NPM);

%% Save width data
save("width_data_NPM.mat",'disp_lower_NPM','disp_upper_NPM','acc_lower_NPM','acc_upper_NPM','vel_lower_NPM',"vel_upper_NPM")