function [U, Udot, Uddot] = Newmark(M, C, K, F, dt, N, Uo, Udoto)
% Newmark Beta Method for Solving Second-Order ODEs
% Solves M*Uddot + C*Udot + K*U = F using Newmark integration.
%
% INPUTS:
% M - Mass matrix (Ndof x Ndof)
% C - Damping matrix (Ndof x Ndof)
% K - Stiffness matrix (Ndof x Ndof)
% F - Force matrix (Ndof x N) [external force at each time step]
% dt - Time step size
% N - Number of time steps
% Uo - Initial displacement vector
% Udoto - Initial velocity vector
%
% OUTPUTS:
% U - Displacement matrix (Ndof x N)
% Udot - Velocity matrix (Ndof x N)
% Uddot - Acceleration matrix (Ndof x N)

[dof, ~] = size(F); % Number of degrees of freedom

% Initialize displacement, velocity, and acceleration matrices
U = zeros(dof, N);
Udot = zeros(dof, N);
Uddot = zeros(dof, N);

% Set initial conditions
U(:, 1) = Uo;
Udot(:, 1) = Udoto;
Uddot(:, 1) = M \ (F(:, 1) - C * Udot(:, 1) - K * U(:, 1));

% Newmark-beta method parameters
delta = 0.6; 
alpha = 0.38;

% Precompute constants for efficiency
a0 = 1 / (alpha * dt^2);
a1 = delta / (alpha * dt);
a2 = 1 / (alpha * dt);
a3 = (1 / (2 * alpha) - 1);
a4 = delta / alpha - 1;
a5 = (dt / 2) * (delta / alpha - 2);
a6 = dt * (1 - delta);
a7 = delta * dt;

% Compute effective system matrix and its inverse
MM = a0 * M + a1 * C + K;
MMM = inv(MM); % Inverse of the effective matrix

% Time integration loop
for i = 1:N-1
    % Compute internal force contributions
    Mt = M * (a0 * U(:, i) + a2 * Udot(:, i) + a3 * Uddot(:, i));
    Ct = C * (a1 * U(:, i) + a4 * Udot(:, i) + a5 * Uddot(:, i));
    
    % Solve for next displacement
    U(:, i+1) = MMM * (F(:, i+1) + Mt + Ct);
    
    % Update velocity and acceleration
    Udot(:, i+1) = Udot(:, i) + dt * ((1 - delta) * Uddot(:, i) ...
        + delta * a0 * (U(:, i+1) - U(:, i) - Udot(:, i) * dt) ...
        - (delta / alpha) * (1/2 - alpha) * Uddot(:, i));
    
    Uddot(:, i+1) = a0 * (U(:, i+1) - U(:, i) - Udot(:, i) * dt) ...
        - (1 / alpha) * (1/2 - alpha) * Uddot(:, i);
end
end