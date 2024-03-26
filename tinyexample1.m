tic
format short
clc;clear;close;
%% === model parameters ===
n = 200; p = 100; n0 = 2;
N = 500;
rho = [0.25 0.5];
mu = zeros(p,1); c = (1:p);
ama = bsxfun(@minus,c,c');
rho = rho(2);
sigma = rho.^(abs(ama));

%% === 'control' parameter ===

a = 100;  % P = 60

%% === True Beta and Phi ===

%% ------------* Case A *------------
Beta = 0.5*[1.0;0;-1.0;0;0;1.0;zeros(p-6,1)];  
PHI =  0.5*[0;1.0;0;-1.0;zeros(p-4,1)]; 

Censorrate = zeros(N,1);
index = find(Beta~=0);
index1 = find(PHI~=0);

%% ----------------------------------------
%     Monte Carlo simulation
% -----------------------------------------------

for iter = 1:N
    iter
    rng(iter)   % Set a random seed to further control reproducibility
    [T,Z,Delta,T1,Z1,Delta1] = survival_data(n,Beta,PHI,mu,sigma,iter); % simulation data

    [ini_beta(:,iter),ini_phi(:,iter)] =...
        initial_beta(n,Z,1e-5,Delta,a);
    
    %MIC_cox(:,iter) = MIC(n,ini_beta(:,iter),Z,1e-5,Delta,a);

end


for iter = 1:N
    iter
    rng(iter)   % Set a random seed to further control reproducibility
    [T,Z,Delta,T1,Z1,Delta1] = survival_data(n,Beta,PHI,mu,sigma,iter);
    Censorrate(iter) = 1-mean(Delta);

    oracle_beta(:,iter) = initial_beta1(n,Z(:,index),1e-5,Delta);
    full_beta(:,iter) = initial_beta1(n,Z,1e-5,Delta);

    [pl(:,iter)] = ista_LQA(n,ini_beta(:,iter),Z,1e-5,Delta,a);  

    [mpl(:,iter),phi(:,iter),seta(:,iter),gama(:,iter),H_T(:,iter),H_C(:,iter)]=...
        MPL(n,n0,mean(ini_beta,2),mean(ini_phi,2),Z1,T1,1e-5,Delta1,a);  


    se_pl(iter) = (pl(:,iter)-Beta)'*sigma*(pl(:,iter)-Beta);    % % MSE
    se_mpl(iter) = (mpl(:,iter)-Beta)'*sigma*(mpl(:,iter)-Beta); % % MSE
    se_phi(iter) = (phi(:,iter)-PHI)'*sigma*(phi(:,iter)-PHI);   % % MSE
    se_oracle(iter)= (oracle_beta(:,iter)-Beta(index))'*sigma(index,index)*(oracle_beta(:,iter)-Beta(index));  % % MSE
    se_full(iter)= (full_beta(:,iter)-Beta)'*sigma*(full_beta(:,iter)-Beta);   % % MSE

    [cov_beta(:,iter),cov_phi(:,iter)] = cov_matrix(mpl(:,iter),phi(:,iter),seta(:,iter)',gama(:,iter)', ...
        T1,Z1,Delta1,n,n0,index,index1);


    PL_cov_beta(:,iter) = PL_cov(n,Z,pl(:,iter),Delta,index);

end

time = toc


%% ---------------------------------------------------------------------------
%      Assessment Criteria :
% ----------------------------------------------------------------------------
% Pcorr: The proportion of all active predictors are selected
% size: The number of selected variables
% MSE: The mean weighted squared error
% F_+: The number of incorrectly selected variables (false negatives)
% F_-: The number of incorrectly excluded variables (false positives)

% -------------------------------------------------------------------------------------
Pcorr_pl = (1/N)*sum((all(pl(index,:))).*(1-any(pl(setdiff(1:1:p, index),:))));
MSE_pl = mean(se_pl);
F_plus_pl = (1/N)*sum(sum(pl(setdiff(1:1:p, index),:)~=0));
F_minus_pl = (1/N)*sum(sum(pl(index,:)==0));
Size_pl = sum(sum(pl(:,:)~=0))/N;

% -------------------------------------------------------------------------------------
Pcorr_mpl = sum((all(mpl(index,:))).*(1-any(mpl(setdiff(1:1:p, index),:))))/N;
MSE_mpl = mean(se_mpl);
F_plus_mpl = sum(sum(mpl(setdiff(1:1:p, index),:)~=0))/N;
F_minus_mpl = sum(sum(mpl(index,:)==0))/N;
Size_mpl = sum(sum(mpl(:,:)~=0))/N;

% -------------------------------------------------------------------------------------
Pcorr_phi = sum((all(phi(index1,:))).*(1-any(phi(setdiff(1:1:p, index1),:))))/N;
MSE_phi = mean(se_phi);
F_plus_phi = sum(sum(phi(setdiff(1:1:p, index1),:)~=0))/N;
F_minus_phi = sum(sum(phi(index1,:)==0))/N;
Size_phi = sum(sum(phi(:,:)~=0))/N;

% -------------------------------------------------------------------------------------
Pcorr_oracle = sum((all(oracle_beta(:,:))))/N;
MSE_oracle = mean(se_oracle);
% F_plus_oracle = sum(sum(oracle_beta(setdiff(index, index),:)~=0))/N;
F_plus_oracle = 0;
F_minus_oracle = sum(sum(oracle_beta(:,:)==0))/N;
Size_oracle = sum(sum(oracle_beta(:,:)~=0))/N;

% -------------------------------------------------------------------------------------
Pcorr_full = sum((all(full_beta(index,:))).*(1-any(full_beta(setdiff(1:1:p, index),:))))/N;
MSE_full = mean(se_full);
F_plus_full = sum(sum(full_beta(setdiff(1:1:p, index),:)~=0))/N;
F_minus_full = sum(sum(full_beta(index,:)==0))/N;
Size_full = sum(sum(full_beta(:,:)~=0))/N;


%% === Output ===
mean(Censorrate)   % Censoring rate

Criteria = [Pcorr_full MSE_full F_plus_full F_minus_full Size_full;
    Pcorr_pl MSE_pl F_plus_pl F_minus_pl Size_pl;
    Pcorr_mpl MSE_mpl F_plus_mpl F_minus_mpl Size_mpl;
    Pcorr_oracle MSE_oracle F_plus_oracle F_minus_oracle length(index)]



%% ===================================================
%                 survival_data()
% ============================================================
function  [T,Z,status,T1,Z1,status1] = survival_data(n,Beta,Phi,mu,sigma,iter)
%% Generating the survival data
p = length(Beta);
ZZ = mvnrnd(mu,sigma,n);  % % n*p covariates
tau = 0.5;  % % Kendall's tau
alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.
u = copularnd('Frank',alph,n);
lamt = 2.3; % tau = 0.5;
% lamt = 2.0;   %tau = 0.2;

Death_time = sqrt(-lamt^2*log(1-u(:,1))./exp(ZZ*Beta));  % % death time
C = -5*log(1-u(:,2))./exp(ZZ*Phi);

statu = (Death_time <= C);
TT = min(Death_time,C);      % % survial time

[T,I] = sort(TT,'descend');  % % sorting the time
% Y = bsxfun(@ge,T,T');      % % at risk process
Z = ZZ(I,:);
status = statu(I);

[T1,I1] = sort(TT,'ascend');  % % sorting the time
% Y = bsxfun(@ge,T1,T1');     % % at risk process
Z1 = ZZ(I1,:);
status1 = statu(I1);

end

%% ===================================================
%                 initial_beta1()
% ============================================================

function initial_beta= initial_beta1(n,Z,r,status)
[~,p] = size(Z);
beta=zeros(p,1);
k = 1; err = 0;  tk = 11;
while k<=1000 && err==0
    % k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;

    beta1 =  beta - L_prime/tk;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;

    k = k+1;
end

beta1 = beta1.*(abs(beta1)>=2*1e-4);

initial_beta = beta1;


end


%% ===================================================
%                 initial_beta()
% ============================================================
function [initial_beta,initial_phi]= initial_beta(n,Z,r,status,a)
[~,p] = size(Z);
beta = zeros(p,1);
phi = zeros(p,1);
opt_BIC = 1e+10;
opt_BIC2 = 1e+10;

n0 = sum(status);
lambda0 = log(n0);


% ============================================================
k = 1; err = 0; 
tk = 18; 

while k<=1000 && err==0
    % k

    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

ini_beta = beta1;

k = 1;err = 0; tk = 18;
beta = ini_beta;
while k<=1000 && err==0
    % k

    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
    u1 = eye(p) + 2*W1/tk;
    beta_tilde = beta - L_prime/tk;
    beta1 = u1\beta_tilde;


    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;

end

beta2 = beta1.*(abs(beta1)>=2.5*1e-4);
opt_beta = beta2;


% ============================================================
k = 1; err1 = 0;
tk = 10; 

while k<=1000 && err1==0
    % k
    L_prime1 = -(sum((1-status).*(Z-cumsum((exp(Z*phi).*Z))./cumsum(exp(Z*phi)))))'/n;
    phi1 =  phi - L_prime1/tk;
    w1 = phi1-phi;
    err1 = norm(w1,2)^2 <= r*norm(phi,2)^2;
    phi = phi1;
    k = k+1;
end
ini_phi = phi1;

k = 1; err1 = 0; tk = 10;
phi =ini_phi;
while k<=1000 && err1==0
    % k
    L_prime1 = -(sum((1-status).*(Z-cumsum((exp(Z*phi).*Z))./ cumsum(exp(Z*phi)))))'/n;

    W2 = diag( (4*a*lambda0*exp(2*a*phi.^2))./((exp(2*a*phi.^2)+1).^2) );
    u2 = eye(p) + 2*W2/tk;
    phi_tilde = phi - L_prime1/tk;
    phi1 = u2\phi_tilde;


    w1 = phi1-phi;
    err1 = norm(w1,2)^2 <= r*norm(phi,2)^2;
    % err1 = sqrt(abs(norm(phi1,2)^2-norm(phi,2)^2)) <= r;
    phi = phi1;
    k = k+1;

end

phi1 = phi1.*(abs(phi1)>=2.5*1e-4);
opt_phi = phi1;


initial_beta = opt_beta;
initial_phi = opt_phi;

end


%% ===================================================
%                 ista_LQA()
% ============================================================
function [opt_beta,opt_theta] = ista_LQA(n,ini_beta,Z,r,status,a)
[~,p] = size(Z);
% beta = ini_beta;
beta = zeros(p,1);
n0 = sum(status);
lambda0 = log(n0);

k = 1; err = 0; tk = 40;
while k<=1000 && err==0
    % k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;
    w = beta1-beta;

    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

ini_beta1 = beta1;

k = 1;err = 0; tk = 4;
beta = ini_beta1;

while k<=1000 && err==0
    % k

    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
    u1 = eye(p) + 2*W1/tk;
    beta_tilde = beta - L_prime/tk;
    beta1 = u1\beta_tilde;


    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;

end

beta2 = beta1.*(abs(beta1)>=2.5*1e-4);
beta = beta2;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1;err=0; tk = 4;
while k<=1000 && err==0
    % k
    W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
    u = eye(p) + 2*W1/tk;

    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    beta_tilde = beta - L_prime/tk;
    beta1 = u\beta_tilde;


    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

beta2 = beta.*(abs(beta)>=5*1e-4);

opt_beta = beta2;

end



%% ===================================================
%                 ista_MIC()
% ============================================================
function opt_beta = MIC(n,ini_beta,Z,r,status,a)
[~,p] = size(Z);
beta = ini_beta;
n0 = sum(status);
lambda0 = log(n0);

k = 1; err = 0; tk = 4;
while k<=1000 && err==0
    % k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;
    w = beta1-beta;
    % err = sqrt(abs(norm(beta1,2)^2-norm(beta,2)^2)) <= r;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

k=1;err=0; tk = 2;
while k<=1000 && err==0
    %   k
    W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
    u = eye(p) + 2/tk*W1;
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    beta_tilde = beta - L_prime/tk;
    beta1 = u\beta_tilde;

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

opt_beta = beta.*(abs(beta)>5*1e-4);

end

%% ===================================================
%                 MPL()
% ============================================================

function [opt_beta,opt_phi,opt_seta,opt_gama,H_T,H_C] = MPL(n,n0,ini_beta,ini_phi,Z,T,r,Delta,a)
[~,p] = size(Z);
m = n/n0;
tau = 0.5;  % % Kendall's tau
alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.

beta = ini_beta;
phi = ini_phi;
lambda0 = log(sum(Delta));  % According to Su et al., 2016
lambda1 = 1e-1*sqrt(n);
lambda2 = 1e-1*sqrt(n);

psix = zeros(n,m);
Psi = zeros(n,m);

%% Discretizing a sample of n survival times into sub-intervals
% with no. of 'binCount' subjects in each sub-interval

%[ID,binwv] = discrBinNA(T,n0);
[ID,~,binwv] = discrBinNA3(T,n0);


%% piecewise constant basis function for the non-parametric baseline hazard
for i=1:m
    psix(((i-1)*n0+1):1:i*n0,i) =1;
end
%% piecewise constant cumulative basis function for the non-parametric baseline hazard
for i = 1:n
    Psi(i,1:ID(i)) = 1*binwv(1:ID(i));
end

%% Initial estimates of theta (piecewise constant estimate of h_{0t}) based on independent censoring assumption
ini_seta = sum( repmat(Delta,1,m).*psix )./( sum(repmat(exp(Z*ini_beta),1,m).*Psi)+1e-6 );

%% Initial estimates of gamma (piecewise constant estimate of h_{0c}) based on independent censoring assumption
ini_gama = sum(repmat(1-Delta,1,m).*psix)./( sum(repmat(exp(Z*ini_phi),1,m).*Psi)+1e-6 );
seta = ini_seta+1e-6;  % 1*m
gama = ini_gama+1e-6;  % 1*m
RT = mat1(psix,T,1e-5); % First order difference -- For piecewise constant

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1; err1=0; err2=0; err3=0; err4=0; 
tk1= 2; tk2= 18;

while k<=1000 && (err1==0 ||err2==0 || err3==0 || err4==0)
    % k
    W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
    u1 = eye(p) + 2*W1/tk1;
    W2 = diag( (4*a*lambda0*exp(2*a*phi.^2))./((exp(2*a*phi.^2)+1).^2) );
    u2 = eye(p) + 2*W2/tk2;


    %%
    H_T = sum( exp(Z*beta).*(seta.*Psi),2 );  % % Cumulative Risk - Failure
    H_C = sum( exp(Z*phi).*(gama.*Psi),2 );   % % Cumulative Risk - Censoring

    S_T = exp(-H_T);  % % Survival Probility - Failure
    S_C = exp(-H_C);  % % Survival Probility - Censoring

    S1 = exp(alph*S_T).*exp(alph*S_C)-exp(alph*S_T)-exp(alph*S_C)+exp(alph);  % % Copula: K(a,b,alph)
    S2 = exp(alph)*exp(alph*S_C)+exp(alph*S_C)-exp(2*alph*S_C)-exp(alph);
    S3 = exp(alph)*exp(alph*S_T)+exp(alph*S_T)-exp(2*alph*S_T)-exp(alph);

    LAM1 = alph*(exp(alph*S_T)./S1).*(Delta.*(S2./(exp(alph*S_T).*(exp(alph*S_C)-1)))+...
        (1-Delta).*(exp(alph)-1)./(exp(alph*S_T)-1));
    LAM2 = (alph./S1).*(Delta.*((exp(alph*S_T)*(exp(alph)-1))./(exp(alph*S_C)-1))+...
        (1-Delta).*(S3./(exp(alph*S_C)-1)));


    DL_beta =  -(sum( (Delta-Delta.*H_T-LAM1.*S_T.*H_T).*Z,1) )'/n;   % % first partial derivative - beta
    DL_phi =   -(sum(((1-Delta)-(1-Delta).*H_C-LAM2.*S_C.*H_C).*Z,1))'/n;  % % first partial derivative - phi
    
 
    beta_tilde = beta - DL_beta/tk1;
    beta1 = u1\beta_tilde;

    w1 = beta1-beta;
    err1 = norm(w1,2)^2 <= r*norm(beta,2)^2;

    beta = beta1;


    phi_tilde = phi - DL_phi/tk2;
    phi1 = u2\phi_tilde;

    w2 = phi1-phi;
    err2 = norm(w2,2)^2 <= r*norm(phi,2)^2;

    phi = phi1;

    %% Armijo's linear search
    uu = 0.6; % uu in (0,1)
    gamm = 0.5; % gamm in (0,0.5)
    sigm = 0.35; % sigm in (0,1)

    t = 1; 

    while (t>0)

        H_T = sum( exp(Z*beta).*(seta.*Psi),2 );  % Cumulative Risk - Failure
        H_C = sum( exp(Z*phi).*(gama.*Psi),2 );   % Cumulative Risk - Censoring


        S_T = exp(-H_T);  % Survival Probility - Failure
        S_C = exp(-H_C);  % Survival Probility - Censoring

        S1 = exp(alph*S_T).*exp(alph*S_C)-exp(alph*S_T)-exp(alph*S_C)+exp(alph);  % Copula: K(a,b,alph)
        S2 = exp(alph)*exp(alph*S_C)+exp(alph*S_C)-exp(2*alph*S_C)-exp(alph);
        S3 = exp(alph)*exp(alph*S_T)+exp(alph*S_T)-exp(2*alph*S_T)-exp(alph);

        %%  Penalty of seta_u and gama_u

        DL_seta =  sum( (Delta.*psix)./(sum(seta.*psix,2))- ...
            (exp(Z*beta).*(Delta+LAM1.*S_T)).*Psi,1)/n - lambda1*seta*RT;   % first partial derivative - seta

        DL_gama =  sum( ((1-Delta).*psix)./(sum(seta.*psix,2))- ...
            (exp(Z*phi).*(1-Delta+LAM2.*S_C)).*Psi,1)/n - lambda2*gama*RT;  % first partial derivative - gama

        LAM11 = alph*(exp(alph*S_T)./S1).*(Delta.*(max(S2,0)./(exp(alph*S_T).*(exp(alph*S_C)-1)))+...
            (1-Delta).*(exp(alph)-1)./(exp(alph*S_T)-1));
        LAM22 = (alph./S1).*(Delta.*((exp(alph*S_T)*(exp(alph)-1))./(exp(alph*S_C)-1))+...
            (1-Delta).*(max(S3,0)./(exp(alph*S_C)-1)));


        Ps1 = sum(( Delta.*exp(Z*beta) ).*Psi + ...
            ((LAM11.*S_T).*exp(Z*beta)).*Psi,1)+lambda1*max(seta*RT,0)+1e-5;
        Ps2 = sum(((1-Delta).*exp(Z*beta)).*Psi + ...
            ((LAM22.*S_C).*exp(Z*phi)).*Psi,1)+lambda2*max(gama*RT,0)+1e-5;

        seta1 = seta + uu*DL_seta*diag(seta./Ps1);
        gama1 = gama + uu*DL_gama*diag(gama./Ps2);


        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S1 = exp(alph*S_T).*exp(alph*S_C)-exp(alph*S_T)-exp(alph*S_C)+exp(alph);
        K1 = (exp(alph*S_T).*(exp(alph*S_C)-1))./S1;
        K2 = (exp(alph*S_C).*(exp(alph*S_T)-1))./S1;
        L_T = log(sum(seta.*psix,2))+Z*beta-H_T+log(K1);
        L_C = log(sum(gama.*psix,2))+Z*phi-H_C+log(K2);

        L_0 = sum(Delta.*L_T+(1-Delta).*L_C)/n;

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H_T = sum( exp(Z*beta).*(seta1.*Psi),2 );  % % Cumulative Risk - Failure
        H_C = sum( exp(Z*phi).*(gama1.*Psi),2 );   % % Cumulative Risk - Censoring


        S_T = exp(-H_T);  % % Survival Probility - Failure
        S_C = exp(-H_C);  % % Survival Probility - Censoring

        S1 = exp(alph*S_T).*exp(alph*S_C)-exp(alph*S_T)-exp(alph*S_C)+exp(alph);
        K1 = (exp(alph*S_T).*(exp(alph*S_C)-1))./S1;
        K2 = (exp(alph*S_C).*(exp(alph*S_T)-1))./S1;
        L_T = log(sum(seta1.*psix,2))+Z*beta-H_T+log(K1);
        L_C = log(sum(gama1.*psix,2))+Z*phi-H_C+log(K2);

        L_new = sum(Delta.*L_T+(1-Delta).*L_C)/n;

      
        if  L_0 + gamm*uu.*[seta1-seta gama1-gama]*([DL_seta DL_gama]') <= L_new
            t = 0;  % % break the loop
            uu_armijo = uu;
        else
            uu = uu*sigm; % % reduce uu, enter the next loop

        end

    end

    seta1 = seta + uu_armijo*DL_seta*diag(seta./Ps1);
    gama1 = gama + uu_armijo*DL_gama*diag(gama./Ps2);

    % seta1 = seta + uu_armijo*(seta1-seta);
    % gama1 = gama + uu_armijo*(gama1-gama);

    err3 = norm(seta1-seta,2)^2 <= 1e-7*norm(seta,2)^2;
    err4 = norm(gama1-gama,2)^2 <= 1e-7*norm(gama,2)^2;

    % err3 = norm(seta1-seta,2) <= 1e-5;
    % err4 = norm(gama1-gama,2) <= 1e-5;

    gama = gama1;
    seta = seta1;

    k = k+1;
end

%%
beta2 = beta1.*(abs(beta1)>=5*1e-4);
phi2 = phi1.*(abs(phi1)>=5*1e-4);

opt_beta = beta2;
opt_phi = phi2;
opt_seta = seta;
opt_gama = gama;


% H_T = sum(opt_seta.*Psi,2);
% H_C = sum(opt_gama.*Psi,2);

end



function [cov_beta, cov_phi] = cov_matrix(beta,phi,seta,gama,T,Z,Delta,n,n0,index,index1)
m = n/n0;
[~,p] = size(Z);
tau = 0.2;  % % Kendall's tau
alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.

psix = zeros(n,m);
Psi = zeros(n,m);
[ID,binwv] = discrBinNA(T,n0);

%% piecewise constant basis function for the non-parametric baseline hazard
for i=1:m
    psix(((i-1)*n0+1):1:i*n0,i) =1;
end
%% piecewise constant cumulative basis function for the non-parametric baseline hazard
for i = 1:n
    Psi(i,1:ID(i)) = 1*binwv(1:ID(i));
end

H_T = sum( exp(Z*beta).*(seta.*Psi),2 );  % % Cumulative Risk - Failure
H_C = sum( exp(Z*phi).*(gama.*Psi),2 );   % % Cumulative Risk - Censoring

S_T = exp(-H_T);  % % Survival Probility - Failure
S_C = exp(-H_C);  % % Survival Probility - Censoring

K1 = zeros(n,1);   K2 = zeros(n,1);    K11 = zeros(n,1);   K211 = zeros(n,1);
K22 = zeros(n,1);  K12 = zeros(n,1);   K21 = zeros(n,1);   K212 = zeros(n,1);
K111 = zeros(n,1); K112 = zeros(n,1);  K221 = zeros(n,1);  K222 = zeros(n,1);


for i=1:n

    K1(i) = ( exp(alph*S_T(i))*(exp(alph*S_C(i))-1) )/( exp(alph*S_T(i))*exp(alph*S_C(i))-...
        exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph) );


    K2(i) = ( exp(alph*S_C(i))*(exp(alph*S_T(i))-1) )/( exp(alph*S_T(i))*exp(alph*S_C(i))-...
        exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph) );


    K11(i) = ( alph*exp(alph*S_T(i))*(exp(alph*S_C(i))*exp(alph)+exp(alph*S_C(i))-exp(2*alph*S_C(i))-exp(alph)) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^2;


    K22(i) = ( alph*exp(alph*S_C(i))*( exp(alph*S_T(i))*exp(alph)+exp(alph*S_T(i))-exp(2*alph*S_T(i))-exp(alph) ) )/...
        ( exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph) )^2;


    K12(i) = ( alph*exp(alph*S_T(i))*exp(alph*S_C(i))*( exp(alph)-1) )/( exp(alph*S_T(i))*exp(alph*S_C(i))-...
        exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph) )^2;


    K21(i) = K12(i);


    K111(i) = ( alph^2*exp(alph*S_T(i))*( exp(alph)*exp(alph*S_C(i))+exp(alph*S_C(i))-exp(2*alph*S_C(i))-exp(alph) )* ...
        (exp(alph*S_T(i))-exp(alph*S_T(i))*exp(alph*S_C(i)))-exp(alph*S_C(i))+exp(alph) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^3;


    K222(i) = ( alph^2*exp(alph*S_C(i))*( exp(alph)*exp(alph*S_T(i))+exp(alph*S_T(i))-exp(2*alph*S_T(i))-exp(alph) )* ...
        (exp(alph*S_C(i))-exp(alph*S_T(i))*exp(alph*S_C(i)))-exp(alph*S_T(i))+exp(alph) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^3;


    K112(i) = ( alph^2*exp(alph*S_T(i))*exp(alph*S_C(i))*( exp(alph)+1-2*exp(alph*S_C(i)) )*...
        ( exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph) )-2*alph^2*exp(alph*S_T(i))*exp(alph*S_C(i))*...
        ( exp(alph*S_C(i))*exp(alph)+exp(alph*S_C(i))-exp(2*alph*S_C(i))-exp(alph) )*( exp(alph*S_T(i))-1 ) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^3;


    K221(i) = ( alph^2*exp(alph*S_T(i))*exp(alph*S_C(i))*( exp(alph)+1-2*exp(alph*S_T(i)) )*...
        ( exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph) )-2*alph^2*exp(alph*S_T(i))*exp(alph*S_C(i))*...
        ( exp(alph*S_T(i))*exp(alph)+exp(alph*S_T(i))-exp(2*alph*S_T(i))-exp(alph) )*( exp(alph*S_C(i))-1 ) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^3;


    K211(i) =  ( alph^2*exp(alph*S_T(i))*exp(alph*S_C(i))*(exp(alph)-1)*( exp(alph*S_T(i))-...
        exp(alph*S_T(i))*exp(alph*S_C(i))- exp(alph*S_C(i)) +exp(alph) ) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^3;

    K212(i) =  ( alph^2*exp(alph*S_T(i))*exp(alph*S_C(i))*(exp(alph)-1)*( exp(alph*S_C(i))-...
        exp(alph*S_T(i))*exp(alph*S_C(i))- exp(alph*S_T(i)) +exp(alph) ) )/...
        (exp(alph*S_T(i))*exp(alph*S_C(i))-exp(alph*S_T(i))-exp(alph*S_C(i))+exp(alph))^3;


end

LAM1 = Delta.*(K11./K1) +(1-Delta).*(K21./K2);

LAM2 = Delta.*(K12./K1) +(1-Delta).*(K22./K2);


L1 = Delta.*(K111.*K1-K11.^2)./(K1.^2)+(1-Delta).*(K112.*K2-K12.^2)./(K2.^2);

L2 = (1-Delta).*(K222.*K2-K22.^2)./(K2.^2)+Delta.*(K221.*K2-K21.^2)./(K1.^2);

L3 = Delta.*( K211.*K1-K11.*K21 )./K1.^2 + (1-Delta).*( (K212.*K2-K12.*K22) ./(K2.^2));

D1 = (Delta-L1.*H_T.*S_T.^2-LAM1.*H_T.*S_T+LAM1.*S_T).*H_T;

D2 = ((1-Delta)-L2.*H_C.*S_C.^2-LAM2.*H_C.*S_C+LAM2.*S_C).*H_C;

D3 = L3.*H_T.*H_C.*S_T.*S_T;

P1 = zeros(p,p);   P2 = zeros(p,p);  P3 = zeros(p,p);

for i=1:p

    for j=1:p

        P1(i,j) = 2*sum(D1.*Z(:,i).*Z(:,j));

        P2(i,j) = 2*sum(D2.*Z(:,i).*Z(:,j));

        P3(i,j) = 2*sum(D3.*Z(:,i).*Z(:,j));

    end

end

tt = [P1 P3; P3 P2];

cov_beta = diag(inv(tt(index,index)));
cov_phi = diag(inv(tt(p+index1,p+index1)));



end



function PL_cov_beta = PL_cov(n,Z,beta,status,index)
%% % * --- ASE:the average of estimated standard error; --- *% % %%
eta = Z*beta;
t = cumsum(exp(eta));  % % n*1
v = cumsum( (exp(eta).*Z) ); % % n*p
for i=1:n
    Cel1(1,i) = {exp(eta(i))*Z(i,:)'*Z(i,:)};
    Cel2(1,i) = {v(i,:)'*v(i,:)};
end
f1 = cumsum(cat(3,Cel1{:}),3);   
f2 = cumsum(cat(3,Cel2{:}),3);  
for i=1:n
    Cel3(1,i) = {f1(:,:,i)};
    Cel4(1,i) = {f2(:,:,i)};
    Cel3(1,i) =  cellfun(@(x) (status(i)/t(i)).*x, Cel3(1, i),'UniformOutput',false);
    Cel4(1,i) =  cellfun(@(x) (status(i)/t(i)^2).*x, Cel4(1,i),'UniformOutput',false);
end
L_primeprime = -2*(sum(cat(3,Cel4{:}),3) - 2*sum(cat(3,Cel3{:}),3));

tt = diag( inv( L_primeprime(index,index) ) );

PL_cov_beta = tt;  % % Cao et al.(2017) Cox-SELO、Zhang and Lu(2007) Adaptive-Cox


end



% % Discretizing a sample of n survival times into sub-intervals
% % with no. of 'binCount' subjects in each sub-interval

%-- binwv: the bin widths
%-- binID: the bin ID of each subject
%-- binedg
function  [binID,binedg,binwv] = discrBinNA3(T, binCount)
[T1,~] = sort(T,'ascend');  % % sorting the time
len = length(T); % sample size
binid = 1:binCount:len;
nbins = length( binid );
binedg = T1(binid);

rank = zeros(size(T1));
for i = 1:len
    rank (i) = find(T1 == T(i), 1);
end

TT = unique(T);
if numel(TT) < numel(T)
    tie =1;
else
    tie=0;
end


if tie == 0

    if binCount == 1
        binwv = binedg(1:nbins) - [0,binedg(1:(nbins-1))];  %binwv: bin widths
        binedg = [0, binedg];
        binID = rank;   %binID: bin IDs
    else
        i = 0;
        while i < nbins-1
            i = i + 1;
            ntied = sum( binedg(i) == binedg );
            if i + ntied  <= nbins
                binedg = [ binedg(1:i);binedg( (i+ntied):nbins)];
            else
                binedg = [ binedg(1:i) ];
            end
            nbins = length( binedg );

        end

        if binedg(nbins) == max(T)
            nbins = nbins - 1;
            binedg = binedg(1:nbins);
        end

        binedg(1) = min(T);
        binedg(nbins + 1) = max(T)+1e-2 *( binedg(nbins)-binedg(nbins-1) );

        binedg(nbins) = binedg(nbins)*( 1-1e-5 );

        binwv = binedg(2:(nbins+1)) - binedg(1:nbins);
        binID = zeros(1,len);
        for i = 1:len
            for j = 1:nbins
                if T(i) < binedg(j+1) && T(i) >= binedg(j)
                    binID(i) = j;
                end
            end
        end
    end

else
    i = 0;
    while i < ( nbins-1 )
        i = i + 1;
        ntied = sum( binedg(i)== binedg );
        if i + ntied  <= nbins
            binedg = [ binedg(1:i);binedg((i+ntied):nbins) ];
        else
            binedg = [ binedg(1:i) ];
        end
        nbins = length( binedg );
    end

    if binedg(nbins) == max(T)
        nbins = nbins - 1;
        binedg = binedg(1:nbins);
    end

    binedg(1) = max(T);
    binedg(nbins+1) = max(T) + 1e-2;
    binedg(nbins) = binedg(nbins)*( 1 -1e-5 );
    binwv = binedg(2:(nbins+1))-binedg(1:nbins);
    binID = zeros(1,len);
    for i = 1:len
        for j = 1:nbins
            if  T(i) < binedg(j+1) && T(i)  >= binedg(j)
                binID(i) = j;
            end
        end
    end

end

end


