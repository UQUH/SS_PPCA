function [K,K_bc,f_vec,Phi,y_EXP,y_HDM_pert,y_ROM,V_global,n,k,...
    y_EXP_no_noise,zeta,num_pertb] = Ex2Model
rng(1042); % Set Random Seed for Reproducibility  

%% Define Parameters
n = 1000; % Number of nodes
n_int = n - 2;  % Interior DOFs (excluding boundary)

% Construct grid points
zeta = linspace(0, 1, n)';  % Equal-distance grid on [0, 1]

% Construct orthogonal matrix S
S = sqrt(2 / (n - 1)) * sin((1:n_int)' * (1:n_int) * pi / (n - 1));

% Construct Lambda matrix (diagonal matrix of eigenvalues)
lambda_k = (2 * pi * (1:n_int)).^2;
lambda_mat = diag(lambda_k);

% Construct Phi matrix
Phi = zeros(n, n_int);
Phi(2:end-1, :) = S;  % Insert S in the interior rows

% Compute stiffness matrix K
K = Phi * lambda_mat * Phi';  % K = Phi * Lambda * Phi^T

% Apply Dirichlet boundary conditions (removing first and last rows/columns)
K_bc = K(2:end-1, 2:end-1);  % Extract submatrix for interior DOFs

%% Define Force Vector  
phi_part = Phi(:, 2:6);  % Extract phi(:, 2) to phi(:, 6)

mu_vec = [0.5, 0.5, 0.5, 0.5, 1]; % Coefficients for force components  

% Compute force vector using basis functions  
g = phi_part * mu_vec';  
% Normalize by maximum absolute value (scale invariant force)
g_max = max(abs(g));
f_vec = g / g_max;

%% Compute Experimental Stiffness matrix  
% Generate random noise vector  
z = randn(1, n_int);  

% Define error percentage  
error_percent = 0.15;  

% Compute scaling factor for error  
c = error_percent * (norm(lambda_k)) / norm(z .* lambda_k);  

% Compute perturbation in eigenvalues  
w = c * z;  
lambda_error = w .* lambda_k;  

lambda_exp = (lambda_k + lambda_error);
% Compute perturbed stiffness matrix  
% K_error = Phi * diag(lambda_error) * Phi';  
% K_true = K + K_error;  
K_true = Phi* diag(lambda_exp) *Phi';
K_true_bc = K_true(2:end-1,2:end-1);

%% Compute Experimental Displacement  
% Solve system for unperturbed experimental displacement
noise_level = 0.05;
y_EXP_no_noise = zeros(n,1);
y_EXP_no_noise(2:end-1) = K_true_bc \ f_vec(2:end-1);  

% Compute standard deviation of displacement  
y_std = std(y_EXP_no_noise);  

%y_error = zeros(n,1);
% Generate random noise for experimental error  
y_error = noise_level * y_std * randn(n, 1);  

% Add noise to displacement to obtain experimental data  
y_EXP = y_EXP_no_noise + y_error;  

%% Compute Perturbed HDM Displacement  
num_pertb = 100; % Number of perturbation realizations  
y_HDM_pert = zeros(n, num_pertb); % Preallocate perturbed displacement matrix  

for i = 1:num_pertb  
    % % Generate random perturbation scaled by error percentage  
    error_HDM = 2.25* error_percent * randn(1, n_int);  

    % Apply perturbation to eigenvalues  
    lambda_error_pert = lambda_k .* (1 + error_HDM);  

    % if any(lambda_error_pert <= 0)
    %     fprintf('lambda is negative')
    % end
    % Construct perturbed stiffness matrix
    K_pert = Phi * diag(lambda_error_pert) * Phi';
    K_pert_bc = K_pert(2:end-1,2:end-1);

    % z = randn(1,n_int);
    % 
    % % Compute scaling factor for error
    % c = 2.3*error_percent * (norm(lambda_k)) / norm(z .* lambda_k);
    % 
    % % Compute perturbation in eigenvalues
    % w = c * z;
    % lambda_error = w .* lambda_k;
    % 
    % lambda_pert = (lambda_k + lambda_error);

    % % Compute perturbed stiffness matrix
    % K_error_pert = Phi * diag(lambda_error) * Phi';
    % % K_pert = K + K_error_pert;
    % K_pert = Phi * diag(lambda_pert) * Phi';
    % K_pert_bc = K_pert(2:end-1,2:end-1);
    % 
    % z = randn(size(mu_vec));
    % f_error = 0.15 * (phi_part * (mu_vec .* z)') / g_max;
    % f_pert = f_vec(2:end-1) + f_pert(2:end-1);
    % Solve system for perturbed HDM displacement  
    y_HDM_pert(2:end-1, i) = K_pert_bc \ f_vec(2:end-1) ;  
end  

%% Compute Reduced-Order Model (ROM) Displacement  
k = 4; % Number of retained modes  

% Perform Singular Value Decomposition (SVD) to obtain reduced basis  
[V_global, ~, ~] = svds(y_HDM_pert, k);  

% Project stiffness matrix onto reduced space  
K_ROM = V_global' * K * V_global;  

% Solve reduced system for generalized coordinates  
q = K_ROM \ (V_global' * f_vec);  

% Compute ROM displacement by projecting back to full space  
y_ROM = V_global * q;  
end