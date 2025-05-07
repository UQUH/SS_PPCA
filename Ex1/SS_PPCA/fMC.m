function f = fMC(beta,V_star,D,K,f_vec,doExp,y_ROM,num_model,n,k,store)
% Persistent variable to store evaluated beta values across function calls
persistent betaEval;
if isempty(betaEval)
    betaEval = [];
end
persistent DTstat
if isempty(DTstat)
    DTstat = table('Size',[0,8], ...
        'VariableNames', {'beta', 'avgLogdosrom', 'sdLogdosrom', 'q25Logdosrom', 'q75Logdosrom', 'medLogdosrom', 'mse','logbeta'},...
        'VariableTypes',{'double','double','double','double','double','double','double','double'});
end
% Set default value for store if not provided
if nargin < 11
    store = true;
end
% If store is true, append the new beta value to betaEval
if store
    betaEval = [betaEval; beta];
end
assignin('base','store_beta_values',betaEval)
% Find the integer part (floor) of beta
beta0 = floor(beta);
w = beta - beta0;  % Calculate the weight based on fractional part

% Find the rows in DTstat for beta0 and beta0 + 1
f0 = DTstat.mse(DTstat.beta == beta0);        % MSE for beta0
f1 = DTstat.mse(DTstat.beta == beta0 + 1);    % MSE for beta0 + 1
y_SROM = zeros(n,num_model);
% Perform linear interpolation
if isempty(f0)
    L2_dist_SROM_ROM = zeros(1000,1);
    for j = 1:1000
        Z = randn(num_model,beta0);
        mat = D.*Z;
        [W,~,~] = svds(mat,k);
        W = V_star*W;
        K_SROM = W'*K*W;
        for kk = 1:100
            q = K_SROM\(W'*f_vec(:,kk));
            % calculate SROM displacement
            y_SROM(:,kk) = W*q;
        end
        L2_dist_SROM_ROM(j,1) =  sqrt(sum((y_SROM - y_ROM).^2,"all"));
    end
    logdosrom = log(L2_dist_SROM_ROM);
    newrow1 = {beta0,mean(logdosrom),std(logdosrom)...
        ,quantile(logdosrom,0.25),quantile(logdosrom,0.75)...
        ,median(logdosrom),mean((L2_dist_SROM_ROM - doExp).^2),log10(beta0)};  % Mean squared error
    DTstat = [DTstat;newrow1];
    % Find the rows in DTstat for beta0
    f0 = DTstat.mse(DTstat.beta == beta0); % MSE for beta0
end
if isempty(f1)
    L2_dist_SROM_ROM = zeros(1000,1);
    for j = 1:1000
        Z = randn(num_model,beta0+1);
        mat = D.*Z;
        [W,~,~] = svds(mat,k);
        W = V_star*W;
        K_SROM = W'*K*W;
        for kk = 1:100
            q = K_SROM\(W'*f_vec(:,kk));
            % calculate SROM displacement
            y_SROM(:,kk) = W*q;
        end
        L2_dist_SROM_ROM(j,1) =  sqrt(sum((y_SROM - y_ROM).^2,"all"));
    end
    logdosrom = log(L2_dist_SROM_ROM);
    newrow2 = {beta0+1,mean(logdosrom),std(logdosrom)...
        ,quantile(logdosrom,0.25),quantile(logdosrom,0.75)...
        ,median(logdosrom),mean((L2_dist_SROM_ROM - doExp).^2),log10(beta0+1)};  % Mean squared error
    DTstat = [DTstat;newrow2];
    f1 = DTstat.mse(DTstat.beta == beta0 + 1);    % MSE for beta0 + 1
end
f = f0 * (1 - w) + f1 * w;
assignin('base','store_DTstat',DTstat)
end