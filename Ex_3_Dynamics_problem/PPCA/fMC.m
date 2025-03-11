function f = fMC(beta,V_r,D_r,K_r,M_r,beta_coff,...
    rank_U,F_r,qoi_ROM,doExp,n,store)
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
if nargin < 13
    store = true;
end
dt = 5e-5;
N = 1537; %number of time step
dof_obs = 23374;
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

% Perform linear interpolation
if isempty(f0)
    L2_dist_SROM_ROM = zeros(1000,1);
    for j = 1:1000
        Z = randn(rank_U,beta0);
        mat = D_r*Z;
        [W,~,~] = svds(mat,n);
        M_SROM = W'*M_r*W;
        K_SROM = W'*K_r*W;
        C_SROM = beta_coff*K_SROM;
        F_SROM = W'*F_r;
        Qo = zeros(n,1);
        Qdoto = zeros(n,1);
        [~,Qdot,~] = Newmark(M_SROM,C_SROM,K_SROM,F_SROM,dt,N,Qo,Qdoto);
        v_SROM = V_r*W*Qdot;
        qoi_SROM = v_SROM(dof_obs,:);
        L2_dist_SROM_ROM(j,1) =  sqrt(sum((qoi_SROM - qoi_ROM).^2,"all"));
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
        Z = randn(rank_U,beta0+1);
        mat = D_r*Z;
        [W,~,~] = svds(mat,n);
        M_SROM = W'*M_r*W;
        K_SROM = W'*K_r*W;
        C_SROM = beta_coff*K_SROM;
        F_SROM = W'*F_r;
        Qo = zeros(n,1);
        Qdoto = zeros(n,1);
        [~,Qdot,~] = Newmark(M_SROM,C_SROM,K_SROM,F_SROM,dt,N,Qo,Qdoto);
        v_SROM = V_r*W*Qdot;
        qoi_SROM = v_SROM(dof_obs,:);
        L2_dist_SROM_ROM(j,1) =  sqrt(sum((qoi_SROM - qoi_ROM).^2,"all"));
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