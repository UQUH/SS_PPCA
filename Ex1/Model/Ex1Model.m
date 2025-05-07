function [K, K_bc, V_global, f_vec, phi_part, y_HDM, y_ROM, n, k, ...
    zeta, num_models] = Ex2Model

rng(1042); % Set Random Seed for Reproducibility

%% Define Parameters
n = 1000;             % Total number of nodes
n_int = n - 2;        % Interior nodes (excluding Dirichlet boundary nodes)

zeta = linspace(0, 1, n)';  % Uniformly spaced grid points in [0, 1]

%% Construct Basis Functions
% Create orthonormal sine basis (excluding boundary)
S = sqrt(2 / (n - 1)) * sin((1:n_int)' * (1:n_int) * pi / (n - 1));

% Create diagonal matrix of eigenvalues (squared frequency terms)
Lambda_mat = diag((2 * pi * (1:n_int)).^2);

% Extend basis to full domain (add zero rows for boundaries)
Phi = zeros(n, n_int); 
Phi(2:end-1, :) = S;

%% Stiffness Matrix Assembly
K = Phi * Lambda_mat * Phi';      % Full stiffness matrix
K_bc = K(2:end-1, 2:end-1);       % Reduced matrix applying Dirichlet BCs

%% Define Parametric Force 
num_models = 100;                 % Number of parameter samples (mu)
mu = lhsdesign(num_models, 5);   % Latin Hypercube Sampling of 5 params

phi_part = Phi(:, 2:6);          % Basis vectors used to define force

% Preallocate memory for force and solution vectors
f_vec = zeros(n, num_models);    % Forcing vectors for each model
y_HDM = zeros(n, num_models);    % High-dimensional model solutions

%% Generate Force and Solve HDM
for i = 1:num_models
    % Generate g_mu by linear combination of basis vectors weighted by mu
    g_mu = phi_part * mu(i, :)';     
    
    % Normalize by maximum absolute value (scale invariant force)
    g_max_mu = max(abs(g_mu));
    f_vec(:, i) = g_mu / g_max_mu;

    % Solve linear system for interior DOFs using reduced stiffness matrix
    y_HDM(2:end-1, i) = K_bc \ f_vec(2:end-1, i);
end

%% Reduced-Order Model (ROM)
k = 4;                                  % Number of reduced basis vectors
[V_global, ~, ~] = svds(y_HDM, k);     % Truncated SVD (efficient for large data)

K_ROM = V_global' * K * V_global;      % Project full stiffness matrix into reduced space
y_ROM = zeros(n, num_models);          % Preallocate ROM solution matrix

for i = 1:num_models
    % Solve reduced system
    q = K_ROM \ (V_global' * f_vec(:, i));
    
    % Map reduced solution back to full space
    y_ROM(:, i) = V_global * q;
end
end