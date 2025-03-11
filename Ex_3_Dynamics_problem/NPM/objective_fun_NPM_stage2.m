%% Stage -2 Optimization objective function for NPM
function J_alpha = objective_fun_NPM_stage2(alpha,M, K, beta_coff,...
    sorted_index, F, qoi_HDM, qoi_ROM, V_global, lambda_v, k_v, L, x)

seed = 42; % Set random seed for reproducibility
rng(seed);

% Define system parameters
gamma = 0.8; % Scaling factor
wj = 0.5; % Weighting factor
dof = 42342; % Number of degrees of freedom
N = 1537; % Number of time steps
dt = 5e-5; % Time step size
k = 10; % Number of modes (or reduced basis size)

%% load up the alpha_opt_stage_1
load("alpha_opt_stage1.mat",'alpha_stage_1');
alpha_vec = [alpha_stage_1(3:end), alpha];

% Create an empty upper triangular matrix
alpha_mat = triu(zeros(k)); 

% Populate the upper triangular matrix (including diagonal) using reshaping
alpha_mat(triu(true(k))) = alpha_vec;

%% Stochastic reduced-order model
v_sim = 100; % Number of simulations
G = zeros(dof, k); % Initialize G matrix
N_0 = 9937; % Number of nodes
G_node = zeros(N_0, 6); % Initialize node-level G matrix
m_dof = 6; % Number of degrees of freedom per node
vp = 20;
v_SROM = zeros(dof,N,v_sim);
for i = 1:v_sim
    for j = 1:k
        for l = 1:m_dof
            a_j_x = pi.*k_v/(L(1)*alpha_stage_1(2));
            a_j_y = pi.*k_v/(L(2)*alpha_stage_1(2));
            a_j_z = pi.*k_v/(L(3)*alpha_stage_1(2));

            a_j_x_i = x(:,1).*a_j_x;
            a_j_y_i = x(:,2).*a_j_y;
            a_j_z_i = x(:,3).*a_j_z;

            G_x = cos(a_j_x_i);
            G_y = cos(a_j_y_i);
            G_z = cos(a_j_z_i);
            S_x = sin(a_j_x_i);
            S_y = sin(a_j_y_i);
            S_z = sin(a_j_z_i);

            C_S_vec = randn(vp,2);
            C_S_mod = sqrt(C_S_vec(:,1).^2 + C_S_vec(:,2).^2);
            C_vec = C_S_vec(:,1)./C_S_mod;
            S_vec = C_S_vec(:,2)./C_S_mod;

            W_bar = sqrt(lambda_v)'.*abs(randn(vp,1));
            W_C = W_bar.*C_vec;
            W_S = W_bar.*S_vec;

            g_X_i = G_x*W_C + S_x*W_S;
            g_Y_i = G_y*W_C + S_y*W_S;
            g_Z_i = G_z*W_C + S_z*W_S;

            G_node(:,l) = (1/sqrt(dof))*g_X_i.*g_Y_i.*g_Z_i;
        end
        G_reshaped = reshape(G_node,[N_0*6,1]);
        G(:,j) = G_reshaped(sorted_index);
    end

    U = G*alpha_mat;
    %A = U; replaced in the code below
    % M = eye(dof);
    D = (V_global'*M*U + U'*M*V_global)/2;
    Z = U - V_global*D;
    H_s_Z = (eye(k) + (alpha_stage_1(1)^2)*Z'*M*Z)^(-0.5);
    W = (V_global + alpha_stage_1(1)*Z)*H_s_Z;
    
    % Compute reduced-order system matrices
    M_SROM = W' * M * W; % Reduced mass matrix
    K_SROM = W' * K * W; % Reduced stiffness matrix
    C_SROM = beta_coff * K_SROM; % Reduced damping matrix
    F_SROM = W' * F; % Reduced force vector

    % Initialize initial conditions
    Qo = zeros(k, 1); % Initial displacement
    Qdoto = zeros(k, 1); % Initial velocity

    % Solve the reduced system using the Newmark-beta method
    [~,Qdot,~] = Newmark(M_SROM,C_SROM,K_SROM,F_SROM,dt,N,Qo,Qdoto);
    
    % Project reduced solutions back to full space
    v_SROM(:,:,i) = W*Qdot;% Velocity
end

% Compute Mean and Standard Deviation for SROM Velocity
mean_v_srom = mean(v_SROM, 3); % Mean velocity across simulations
std_v_srom = std(v_SROM, 0, 3); % Standard deviation of velocity across simulations

% Compute standard deviation for ROM velocity scaled by gamma
std_v_ROM = gamma * abs(qoi_HDM - qoi_ROM);

% Compute Error Metrics
J_mean = sum((qoi_HDM - mean_v_srom).^2, "all") / sum(qoi_HDM.^2, "all"); % Mean error metric
J_std = sum((std_v_ROM - std_v_srom).^2, "all") / sum(std_v_ROM.^2, "all"); % Standard deviation error metric

% Compute Final Cost Function
J_alpha = wj * J_mean + (1 - wj) * J_std; % Weighted sum of mean and std error

end