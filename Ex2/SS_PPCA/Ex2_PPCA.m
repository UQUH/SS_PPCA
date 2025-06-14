%% Linear static problem: SS-PPCA
rng(42); % Set random seed for reproducibility

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
[K, K_bc, f_vec, phi, y_EXP, y_HDM, y_ROM, V_global, n, k, ...
    y_true_no_noise, x,num_pertb] = Ex2Model;

% Select observation points
points = 1:50:1000; 
points(1) = [];
x_observe = x(points); % Extract observed x-coordinates  
% Extract observed values at selected points
y_EXP_obs = y_EXP(points); % Observed ground truth values
y_ROM_obs = y_ROM(points); % Observed ROM values

%% Perfrom SVD on perturbed HDM Models data
[V_star,D,~] = svd(y_HDM,"econ"); %economy-sized SVD
D = diag(D);% transform into column vector

%% Optimization for Hyperparameter
% Compute L2 distance between EXP & ROM Observed values
dist_EXP_ROM = sqrt(sum((y_EXP_obs-y_ROM_obs).^2,"all")); 

tic; % Start timing
boundBeta = [k, 20]; % Define search bounds for beta optimization

options = optimset('TolX', 0.01); % Set optimization options (tolerance for fminbnd)

% Optimize beta using fminbnd
[betaOpt, ~] = fminbnd(@(beta) fMC(beta,V_star,D,K,f_vec,...
    dist_EXP_ROM,y_ROM_obs,num_pertb,points,k), ...
    boundBeta(1), boundBeta(2), options);

beta_opt = round(betaOpt);    % Corresponding optimal beta

% Compute and display execution time
Opt_time_PPCA = toc; % End timing
fprintf('Hyper paramaeter optimization completed in %.2f seconds.\n', Opt_time_PPCA);

%% Compute Unperturbed HDM Displacement  
Disp_HDM = zeros(n,1);
Disp_HDM(2:end-1,1) = K_bc \ f_vec(2:end-1); % Solve for HDM displacement using the system matrix  

%% Custom Plot Properties  
set(0, 'DefaultAxesFontSize', 12); % Set default font size for axes  
set(0, 'DefaultAxesFontName', 'Times'); % Set default font to Times  
set(0, 'DefaultAxesTickLabelInterpreter', 'latex'); % Use LaTeX formatting for tick labels  
set(0, 'DefaultLegendInterpreter', 'latex'); % Use LaTeX formatting for legends  
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % Ensure global setting for LaTeX tick labels  

% Define figure dimensions  
width_plot = 540; % Figure width in points  
height_plot = width_plot / 3; % Maintain aspect ratio for height  

% Define plot line thickness  
th = 1.4; % Line thickness  

%% To store figures
% Create the folder if it doesn't exist
if ~exist(fullfile(pwd, 'Results'), 'dir')
    mkdir(fullfile(pwd, 'Results'));
end

%% Plot: Experimental Displacement vs HDM  
figure('Units', 'points', 'Position', [0, 0, width_plot, height_plot]); % Create figure  
p0 = plot(x, y_true_no_noise, 'b', 'LineWidth', th);  % Plot experimental displacement 
hold on;  
% Scatter plot for observed experimental data  
p1 = scatter(x_observe, y_EXP_obs, 'b', 'LineWidth', th);  

p2 = plot(x, Disp_HDM, 'k', 'LineWidth', th);  % Plot HDM displacement 
%p2.Color = '#ee3a2b'; % Red color for HDM  
xlabel('X');  
ylabel('Displacement');  
hold off
ylim([-1.2 2.2]*1e-3);   
legend([p1, p0, p2],'Ground truth (observed)','Ground truth (no noise)', 'HDM', 'Box', 'off', 'Location', 'northeast');  

filename = "EXP_HDM_Ex2.pdf";  
output_path = fullfile(pwd, 'Results', filename);% Combine into a full relative path
exportgraphics(gcf, output_path, 'ContentType', 'vector');% Export the figure

%% Compute Stochastic Reduced-Order Basis and Solve for SROM Response  
num_sim = 1e4; % Number of Monte Carlo simulations  
y_SROM = zeros(n, num_sim); % Preallocate SROM displacement matrix  

for i = 1:num_sim
    Z = randn(num_pertb, beta_opt);  % Generate random matrix
    mat = D.* Z;
    [W, ~, ~] = svds(mat, k);  % Perform truncated SVD to get reduced basis 
    W = V_star*W;
    K_SROM = W' * K * W;  % Project system matrix onto reduced space
    q = K_SROM \ (W' * f_vec);  % Solve reduced system
    y_SROM(:, i) = W * q;  % Project back to full space
end  

% Compute mean response of SROM  
mean_y_srom = mean(y_SROM, 2);  

%% Plot: SROM Prediction with 95% Prediction Interval  
pc = 0.95; % 95% confidence level  
alpha = 1 - pc;  

% Compute lower and upper quantiles for the SROM displacement  
y_lower_PPCA = quantile(y_SROM', alpha / 2);  
y_upper_PPCA = quantile(y_SROM', 1 - alpha / 2);  

% Create figure  
figure('Units', 'points', 'Position', [0, 0, width_plot, height_plot]);  

% Plot 95% prediction interval as a shaded region  
h = fill([x', flip(x')], [y_lower_PPCA, flip(y_upper_PPCA)], 'c', 'FaceAlpha', 1);  
h.FaceColor = '#a6cce3'; % Light blue color for the interval  
h.EdgeColor = "none";  

hold on;  

% Scatter plot for observed experimental data  
p1 = scatter(x_observe, y_EXP_obs, 'b', 'LineWidth', th);  

% Plot ground truth without noise  
p0 = plot(x, y_true_no_noise, 'b', 'LineWidth', th);  

% Plot HDM displacement with custom color  
p2 = plot(x, Disp_HDM, 'k', 'LineWidth', th);
%p2.Color = '#ee3a2b'; % Red color for HDM  

% Plot mean SROM displacement  
p3 = plot(x, mean_y_srom, 'LineWidth', th);
p3.Color = "#1e78b3"; % Blue color for SROM mean  

% Set axis labels and ticks  
xlabel('X');  
xticks(0:0.2:1);  
ylabel('Displacement');  
ylim([-1.3 2.25]*1e-3)
hold off

% Add legend with LaTeX formatting  
legend([p1, p0, p2, p3, h], {'Ground truth (observed)', 'Ground truth (no noise)', ...
    'HDM', 'SROM mean', 'SROM $95\%$ PI'}, 'Location', 'northeast','Box', 'off');  

filename = "Prediction_PPCA_Ex2.pdf";  
output_path = fullfile(pwd, 'Results', filename);% Combine into a full relative path
exportgraphics(gcf, output_path, 'ContentType', 'vector');% Export the figure

%% Calculate the coverage and width for displacment
coverage_PPCA = coverage(y_upper_PPCA',y_lower_PPCA',y_true_no_noise);
width_PPCA = mean(y_upper_PPCA - y_lower_PPCA);

%% Calculate the coverage for observed displacment
y_upper_observed = y_upper_PPCA(points);
y_lower_observed = y_lower_PPCA(points);
coverage_observed_data = coverage(y_upper_observed',y_lower_observed',y_EXP_obs);

%% Using real value

%% If need to get the real value of beta 
tic;
% if betaOpt > beta_opt
%     boundBeta = [beta_opt, beta_opt+2];
% else
%     boundBeta = [beta_opt-2,beta_opt];
% % end
clear fmc
boundBeta = [4,6];
options = optimset('TolX', 1e-16); % Set optimization options (tolerance for fminbnd)

% Optimize beta using fminbnd
[betaOpt_real, fval] = fminbnd(@(beta) fMC_real(beta,V_star,D,K,f_vec,...
    dist_EXP_ROM,y_ROM_obs,num_pertb,points,k), ...
    boundBeta(1), boundBeta(2), options);

% Compute and display execution time
Opt_time_PPCA_real = toc; % End timing
fprintf('Hyper paramaeter optimization completed in %.2f seconds.\n', Opt_time_PPCA_real);

%% Compute Stochastic Reduced-Order Basis and Solve for SROM Response  
y_SROM_real = zeros(n, num_sim); % Preallocate SROM displacement matrix  
beta_opt_whole = floor(betaOpt_real);
beta_opt_frac = betaOpt_real - beta_opt_whole;
for i = 1:num_sim
    Z = randn(num_pertb, beta_opt_whole+1);  % Generate random matrix
    Z(:,end) = Z(:,end)*beta_opt_frac;
    mat = D.* Z;
    [W, ~, ~] = svds(mat, k);  % Perform truncated SVD to get reduced basis 
    W = V_star*W;
    K_SROM = W' * K * W;  % Project system matrix onto reduced space
    q = K_SROM \ (W' * f_vec);  % Solve reduced system
    y_SROM_real(:, i) = W * q;  % Project back to full space
end  

% Compute mean response of SROM  
mean_y_srom_real = mean(y_SROM_real, 2);  

%% Plot: SROM Prediction with 95% Prediction Interval  
pc = 0.95; % 95% confidence level  
alpha = 1 - pc;  

% Compute lower and upper quantiles for the SROM displacement  
y_lower_PPCA_real = quantile(y_SROM_real', alpha / 2);  
y_upper_PPCA_real = quantile(y_SROM_real', 1 - alpha / 2);  

% Create figure  
f = figure('Units', 'points', 'Position', [0, 0, width_plot, height_plot]);  

% Plot 95% prediction interval as a shaded region  
h = fill([x', flip(x')], [y_lower_PPCA_real, flip(y_upper_PPCA_real)], 'c', 'FaceAlpha', 1);  
h.FaceColor = '#a6cce3'; % Light blue color for the interval  
h.EdgeColor = "none";  

hold on;  

% Scatter plot for observed experimental data  
p1 = scatter(x_observe, y_EXP_obs, 'b', 'LineWidth', th);  

% Plot ground truth without noise  
p0 = plot(x, y_true_no_noise, 'b', 'LineWidth', th);  

% Plot HDM displacement with custom color  
p2 = plot(x, Disp_HDM, 'k', 'LineWidth', th);
%p2.Color = '#ee3a2b'; % Red color for HDM  

% Plot mean SROM displacement  
p3 = plot(x, mean_y_srom_real, 'LineWidth', th);
p3.Color = "#1e78b3"; % Blue color for SROM mean  

% Set axis labels and ticks  
xlabel('X');  
xticks(0:0.2:1);  
ylabel('Displacement');  
ylim([-1.3 2.25]*1e-3)
hold off

% Add legend with LaTeX formatting  
legend([p1, p0, p2, p3, h], {'Ground truth (observed)', 'Ground truth (no noise)', ...
    'HDM', 'SROM mean', 'SROM $95\%$ PI'}, 'Location', 'northeast','Box', 'off');  

filename = "Prediction_PPCA_Ex2_real.pdf";  
output_path = fullfile(pwd, 'Results', filename);% Combine into a full relative path
exportgraphics(gcf, output_path, 'ContentType', 'vector');% Export the figure

%% Calculate the coverage and width for displacment
coverage_PPCA_real = coverage(y_upper_PPCA_real',y_lower_PPCA_real',y_true_no_noise);
width_PPCA_real = mean(y_upper_PPCA_real - y_lower_PPCA_real);

%% Calculate the coverage for observed displacment
y_upper_observed_real = y_upper_PPCA_real(points);
y_lower_observed_real = y_lower_PPCA_real(points);
coverage_observed_data_real= coverage(y_upper_observed_real',y_lower_observed_real',y_EXP_obs);
