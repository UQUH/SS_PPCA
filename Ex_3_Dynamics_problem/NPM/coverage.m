function coverage_percent = coverage(y_upper, y_lower, qoi_HDM)
% Computes the coverage percentage of HDM quantities within the prediction bounds.
%
% INPUTS:
% y_upper - Upper prediction bound (vector)
% y_lower - Lower prediction bound (vector)
% qoi_HDM - High-Dimensional Model (HDM) quantities of interest (vector)
%
% OUTPUT:
% coverage_percent - Percentage of time steps where qoi_HDM lies within [y_lower, y_upper]

% Define the number of initial time steps to exclude
subset = 200; % Remove the first 200 time steps
N = length(qoi_HDM); % Total number of time steps

% Extract subset after the removed initial steps
qoi_subset = qoi_HDM(subset+1:end);
y_upper_subset = y_upper(subset+1:end);
y_lower_subset = y_lower(subset+1:end);

% Compute coverage violations
upper_violations = sum(qoi_subset > y_upper_subset); % Count exceedances above upper bound
lower_violations = sum(qoi_subset < y_lower_subset); % Count values below lower bound

% Compute coverage percentage
coverage_percent = 100 * (1 - (upper_violations + lower_violations) / (N - subset));
end