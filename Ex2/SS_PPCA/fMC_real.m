function f = fMC_real(beta,V_star,D,K,f_vec,doExp,y_ROM_obs,num_pertb,points,k,store)
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
beta_whole = floor(beta);
beta_frac = beta - beta_whole;

% Find the rows in DTstat for beta
f = DTstat.mse(DTstat.beta == beta);        % MSE for beta

if isempty(f)
    L2_dist_SROM_ROM = zeros(1e5,1);
    for j = 1:1e5
        Z = randn(num_pertb,beta_whole+1);
        Z(:,end) = Z(:,end)*beta_frac;
        mat = D.*Z;
        [W,~,~] = svds(mat,k);
        W = V_star*W;
        K_SROM = W'*K*W;
        q = K_SROM\(W'*f_vec);
        % calculate SROM displacement
        y_SROM = W*q;
        y_SROM_obs = y_SROM(points);
        L2_dist_SROM_ROM(j,1) =  sqrt(sum((y_SROM_obs - y_ROM_obs).^2,"all"));
    end
    logdosrom = log(L2_dist_SROM_ROM);
    newrow1 = {beta,mean(logdosrom),std(logdosrom)...
        ,quantile(logdosrom,0.25),quantile(logdosrom,0.75)...
        ,median(logdosrom),mean((L2_dist_SROM_ROM - doExp).^2),log10(beta)};  % Mean squared error
    DTstat = [DTstat;newrow1];
    % Find the rows in DTstat for beta0
    f = DTstat.mse(DTstat.beta == beta); % MSE for beta0
end
assignin('base','store_DTstat',DTstat)
end