%% Generate HDM and ROM Response Using System Matrices

%% Custom Plot Specifications
width_plot = 1400;
height_plot = 700;

% Set default fonts and interpreters
set(0, 'DefaultTextFontSize', 12); 
set(0, 'DefaultAxesFontSize', 12); 
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Sans Serif');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Load External Force Data
load("force.mat");
N = 1537; % Number of time steps
dt = 5e-5; % Time step size
t = 1000 * dt * (0:N-1); % Time vector (ms)

%% Load Mass, Stiffness, and Damping Matrices for HDM
beta_coff = 6.366e-06;
load("system_matrices.mat");

%% Solve HDM Using Newmark Method
tic;
[dof, ~] = size(K); % Get degrees of freedom
Uo = zeros(dof, 1); % Initial displacement
Udoto = zeros(dof, 1); % Initial velocity

% Solve full-order system
[U, Udot, Uddot] = Newmark(M, C, K, F, dt, N, Uo, Udoto);
save("displacement.mat", "Uddot", "Udot", "U");
toc;

%% Reduced-Order Model (ROM) Computation
tic;
[V, ~, ~] = svd(U, "econ"); % Perform economy-sized SVD
k = 10; % Number of retained modes
V_global = V(:, 1:k); % Reduced-order basis

% Project system matrices onto reduced space
M_ROM = V_global' * M * V_global;
K_ROM = V_global' * K * V_global;
C_ROM = beta_coff * K_ROM;
F_ROM = V_global' * F;

% Solve reduced system
qo = zeros(k, 1); % Initial displacement in reduced space
qdoto = zeros(k, 1); % Initial velocity in reduced space
[q, qdot, qddot] = Newmark(M_ROM, C_ROM, K_ROM, F_ROM, dt, N, qo, qdoto);

% Project back to full space
y_ROM = V_global * q;
v_ROM = V_global * qdot;
a_ROM = V_global * qddot;

toc;
save("displacement_rom_k_10.mat", "a_ROM", "v_ROM", "y_ROM");

%% Plot HDM vs ROM Velocity at Observed DOF
dof_obs = 23374;
th = 2; % Line thickness

f = figure('Color', [1 1 1], 'units', 'points', 'position', [0, 0, width_plot, height_plot]); 
hold on;
plot(t, Udot(dof_obs, :), 'r', 'LineWidth', th);
plot(t, v_ROM(dof_obs, :), 'k', 'LineWidth', th);
hold off;

xlabel('Time (ms)');
ylabel('Velocity in X (in/s)');
legend('HDM', 'ROM');
exportgraphics(gcf, 'HDM_ROM.pdf', 'ContentType', 'vector');

%% Specify Quantity of Interest (QoI)
store_node_states = zeros(18, N);
store_node_states(1:6, :) = U(dof_obs:dof_obs+5, :); % Displacement
store_node_states(7:12, :) = Udot(dof_obs:dof_obs+5, :); % Velocity
store_node_states(13:18, :) = Uddot(dof_obs:dof_obs+5, :); % Acceleration

store_node_states_ROM = zeros(18, N);
store_node_states_ROM(1:6, :) = y_ROM(dof_obs:dof_obs+5, :); % Displacement
store_node_states_ROM(7:12, :) = v_ROM(dof_obs:dof_obs+5, :); % Velocity
store_node_states_ROM(13:18, :) = a_ROM(dof_obs:dof_obs+5, :); % Acceleration

% Save QoI data
save('store_node_states_ROM.mat', 'store_node_states_ROM');
save('store_node_states_HDM.mat', 'store_node_states');

% Clear unnecessary variables
clear U Udot Uddot y_ROM v_ROM a_ROM;