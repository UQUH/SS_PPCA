function coverage_percent = coverage(y_upper, y_lower, y_EXP)
% Computes the coverage percentage of  EXP quantities
% within the prediction bounds.
%
% INPUTS:
% y_upper - Upper prediction bound (vector)
% y_lower - Lower prediction bound (vector)
% y_EXP - EXP displacement
%
% OUTPUT:
% coverage_percent - Percentage of times where y_EXP lies
% within [y_lower, y_upper]

N = length(y_EXP); % Total number of time steps

% Compute coverage violations
upper_violations = sum(y_EXP > y_upper); % Count exceedances above upper bound
lower_violations = sum(y_EXP < y_lower); % Count values below lower bound

% Compute coverage percentage
coverage_percent = 100 * (1 - (upper_violations + lower_violations) / (N));
end