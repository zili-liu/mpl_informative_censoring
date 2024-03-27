tic
format short
clc;clear;close;
%% === model parameters ===
n = 200; p = 20; n0 = 2;
N = 500;
rho = [0.25 0.5];
mu = zeros(p,1); c = (1:p);
ama = bsxfun(@minus,c,c');
rho = rho(2);
sigma = rho.^(abs(ama));

%% === 'control' parameter ===
a = 120;  
tau = 0.5;

%% === True parameters ===
%% ------------* Case A *------------
Beta = 0.5*[1.0;0;-1.0;0;0;1.0;zeros(p-6,1)];  
PHI =  0.5*[0;1.0;0;-1.0;zeros(p-4,1)];  

Censorrate = zeros(N,1);
index = find(Beta~=0);
index1 = find(PHI~=0);

%% ----------------------------------------
%     Monte Carlo simulation
% -----------------------------------------------

for iter = 1:200
    iter
    rng(iter)   % % Set a random seed
    [T,Z,Delta,T1,Z1,Delta1] = survival_data(n,Beta,PHI,mu,sigma,tau,iter);
    Censorrate1(iter) = 1-mean(Delta);

    [ini_beta(:,iter),ini_phi(:,iter)] =...
        initial_beta(n,Z,1e-5/2,Delta,a);
    % MIC_cox(:,iter) = MIC(n,ini_beta(:,iter),Z,1e-5,Delta,a);
end
mean(Censorrate1)
beta0 = mean(ini_beta,2);
phi0 = mean(ini_phi,2);

for iter = 1:N
    iter
    % rng(iter)   % Set a random seed to further control reproducibility
    [T,Z,Delta,T1,Z1,Delta1] = survival_data(n,Beta,PHI,mu,sigma,tau,iter);
    Censorrate(iter) = 1-mean(Delta);

     %% === save simulation data ===
    save('D:\Informative censoring\simulation revision\Spline\mat\Delta.mat','Delta1','-append');
    save('D:\Informative censoring\simulation revision\Spline\mat\T.mat','T1','-append');
    save('D:\Informative censoring\simulation revision\Spline\mat\Z.mat','Z1','-append');
    save('D:\Informative censoring\simulation revision\Spline\mat\Beta0.mat','beta0','-append');
    save('D:\Informative censoring\simulation revision\Spline\mat\Phi0.mat','phi0','-append');

    %% === calling the R ===
    Rpath = 'D:\Software\R-4.2.0\bin';
    RscriptFileName = 'D:\Informative censoring\simulation revision\Spline\Spline_MPL.R';
    RunRcode(RscriptFileName, Rpath);

   

    %% === read csv data ===
    thet = csvread('D:\Informative censoring\simulation revision\Spline\theta_ini.csv', 1, 1);
    gamm = csvread('D:\Informative censoring\simulation revision\Spline\gamma_ini.csv', 1, 1);
    RT = csvread('D:\Informative censoring\simulation revision\Spline\RT.csv', 1, 1);
    Psi = csvread('D:\Informative censoring\simulation revision\Spline\Psi.csv', 1, 1);
    psix = csvread('D:\Informative censoring\simulation revision\Spline\Psix.csv', 1, 1);


    [mpl(:,iter),phi(:,iter),seta(:,iter),gama(:,iter),H_T(:,iter),H_C(:,iter)]=...
        MPL(n,n0,beta0,phi0,Z1,T1,1e-5,Delta1,a,tau);  

    [mpl1(:,iter),phi1(:,iter),seta1(:,iter),gama1(:,iter),H_T1(:,iter),H_C1(:,iter)] =...
        MPL_mspline(n,n0,beta0,phi0,Z1,1e-5,psix,Psi,thet',gamm',RT,Delta1,a,tau);

    se_mpl(iter) = (mpl(:,iter)-Beta)'*sigma*(mpl(:,iter)-Beta);    % % MSE
    se_mpl1(iter) = (mpl1(:,iter)-Beta)'*sigma*(mpl1(:,iter)-Beta); % % MSE


end


% -------------------------------------------------------------------------------------

Pcorr_mpl = sum((all(mpl(index,:))).*(1-any(mpl(setdiff(1:1:p, index),:))))/N;
MSE_mpl = mean(se_mpl);
F_plus_mpl = sum(sum(mpl(setdiff(1:1:p, index),:)~=0))/N;
F_minus_mpl = sum(sum(mpl(index,:)==0))/N;
Size_mpl = sum(sum(mpl(:,:)~=0))/N;


% -------------------------------------------------------------------------------------

Pcorr_mpl1 = sum((all(mpl1(index,:))).*(1-any(mpl1(setdiff(1:1:p, index),:))))/N;
MSE_mpl1 = mean(se_mpl1);
F_plus_mpl1 = sum(sum(mpl1(setdiff(1:1:p, index),:)~=0))/N;
F_minus_mpl1 = sum(sum(mpl1(index,:)==0))/N;
Size_mpl1 = sum(sum(mpl1(:,:)~=0))/N;


%% output
Criteria = [
    Pcorr_mpl MSE_mpl F_plus_mpl F_minus_mpl Size_mpl;
    Pcorr_mpl1 MSE_mpl1 F_plus_mpl1 F_minus_mpl1 Size_mpl1]

mean(Censorrate)   %% Censoring rate

[mean(mpl(index,:),2),mean(mpl1(index,:),2)]


%% ===================================================
%                 survival_data()
% ============================================================
function  [T,Z,status,T1,Z1,status1] = survival_data(n,Beta,Phi,mu,sigma,tau,iter)
%% Generating the survival data
p = length(Beta);
ZZ = mvnrnd(mu,sigma,n);  % % n*p covarites
%tau = 0.2;  % % Kendall's tau
alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.
u = copularnd('Frank',alph,n);
lamt = 2.3; % tau = 0.5;

Death_time = sqrt(-lamt^2*log(1-u(:,1))./exp(ZZ*Beta));  % % death time
C = -5*log(1-u(:,2))./exp(ZZ*Phi);

statu = (Death_time <= C);
TT = min(Death_time,C);     % % survial time

[T,I] = sort(TT,'descend');  % % sorting the time
% Y = bsxfun(@ge,T,T');     % % at risk process
Z = ZZ(I,:);
status = statu(I);

[T1,I1] = sort(TT,'ascend');  % % sorting the time
% Y = bsxfun(@ge,T1,T1');     % % at risk process
Z1 = ZZ(I1,:);
status1 = statu(I1);

end


function  [T,Z,status,T1,Z1,status1] = survival_data1(n,Beta,Phi,mu,sigma,tau,iter)
%% Generating the survival data
p = length(Beta);
ZZ = mvnrnd(mu,sigma,n);  % % n*p covariates
%tau = 0.2;  % % Kendall's tau
alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.
u = copularnd('Frank',alph,n);
% lamt = 2.3; % tau = 0.5;
% lamt = 2.0;   %tau = 0.2;
lamt = 2.1;   %tau = 0.2;
Death_time = sqrt(-lamt^2*log(1-u(:,1))./exp(ZZ*Beta));  % % death time
C = -5*log(1-u(:,2))./exp(ZZ*Phi);

statu = (Death_time <= C);
TT = min(Death_time,C);     % % survial time

[T,I] = sort(TT,'descend');  % % sorting the time
% Y = bsxfun(@ge,T,T');     % % at risk process
Z = ZZ(I,:);
status = statu(I);

[T1,I1] = sort(TT,'ascend');  % % sorting the time
% Y = bsxfun(@ge,T1,T1');     % % at risk process
Z1 = ZZ(I1,:);
status1 = statu(I1);

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
k = 1; err = 0; tk = 60;
while k<=1000 && err==0
    %     k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

ini_beta = beta1;

k = 1;err = 0; tk = 4;
beta = ini_beta;
while k<=1000 && err==0
    %k

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

phi1 = phi1.*(abs(phi1)>=2*1e-4);
opt_phi = phi1;

initial_beta = opt_beta;
initial_phi = opt_phi;

end


%% ===================================================
%                 MPL()
% ============================================================

function [opt_beta,opt_phi,opt_seta,opt_gama,H_T,H_C] = MPL(n,n0,ini_beta,ini_phi,Z,T,r,Delta,a,tau)
[~,p] = size(Z);
m = n/n0;
% tau = 0.2;  % % Kendall's tau
alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.

beta = ini_beta;
phi = ini_phi;
lambda0 = log(sum(Delta));  % According to Su et al., 2016

% lambda1 = 1e-1*sqrt(n);
% lambda2 = 1e-1*sqrt(n);
lambda1 = 1e+3*sqrt(n);
lambda2 = 1e+3*sqrt(n);

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

%% Initial estimates of theta (piecewise constant estimate of h_{0t}) based on independent censoring assumption
ini_seta = sum( repmat(Delta,1,m).*psix )./( sum(repmat(exp(Z*ini_beta),1,m).*Psi)+1e-6 );

%% Initial estimates of gamma (piecewise constant estimate of h_{0c}) based on independent censoring assumption
ini_gama = sum(repmat(1-Delta,1,m).*psix)./( sum(repmat(exp(Z*ini_phi),1,m).*Psi)+1e-6 );

seta = ini_seta+1e-6;  % 1*m
gama = ini_gama+1e-6;  % 1*m
RT = mat1(psix,T,1e-5); % First order difference -- For piecewise constant

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1; err1=0; err2=0; err3=0; err4=0; 
tk1= 46; tk2= 10;

while k<=1000 && (err1==0 ||err2==0 || err3==0 || err4==0)
    %k
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

    DL_beta =  -(sum( (Delta-Delta.*H_T-LAM1.*S_T.*H_T).*Z,1) )'/n;  % % first partial derivative - beta
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
    uu = 0.6; % % taking value between 0 and 1
    gamm = 0.5; % % taking value between 0 and 0.5
    sigm = 0.35; % % taking value between 0 and 1

    t = 1; % % marking the loop

    while (t>0)

        H_T = sum( exp(Z*beta).*(seta.*Psi),2 );  % % Cumulative Risk - Failure
        H_C = sum( exp(Z*phi).*(gama.*Psi),2 );   % % Cumulative Risk - Censoring

        S_T = exp(-H_T);  % % Survival Probility - Failure
        S_C = exp(-H_C);  % % Survival Probility - Censoring

        S1 = exp(alph*S_T).*exp(alph*S_C)-exp(alph*S_T)-exp(alph*S_C)+exp(alph);  % % Copula: K(a,b,alph)
        S2 = exp(alph)*exp(alph*S_C)+exp(alph*S_C)-exp(2*alph*S_C)-exp(alph);
        S3 = exp(alph)*exp(alph*S_T)+exp(alph*S_T)-exp(2*alph*S_T)-exp(alph);

        DL_seta =  sum( (Delta.*psix)./(sum(seta.*psix,2))- ...
            (exp(Z*beta).*(Delta+LAM1.*S_T)).*Psi,1)/n - lambda1*seta*RT;   % % first partial derivative - seta

        DL_gama =  sum( ((1-Delta).*psix)./(sum(seta.*psix,2))- ...
            (exp(Z*phi).*(1-Delta+LAM2.*S_C)).*Psi,1)/n - lambda2*gama*RT;  % % first partial derivative - gama

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
            t = 0; % % break the loop
            uu_armijo = uu;
        else
            uu = uu*sigm; % % reducing uu, enter the next loop
        end
    end

    seta1 = seta + uu_armijo*DL_seta*diag(seta./Ps1);
    gama1 = gama + uu_armijo*DL_gama*diag(gama./Ps2);
    err3 = norm(seta1-seta,2)^2 <= 1e-7*norm(seta,2)^2;
    err4 = norm(gama1-gama,2)^2 <= 1e-7*norm(gama,2)^2;
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

end



%% ===================================================
%                 MPL_mspline()
% ============================================================

function [opt_beta,opt_phi,opt_seta,opt_gama,H_T,H_C] = MPL_mspline(n,n0,ini_beta,ini_phi,Z,r,psix,Psi,seta,gama,RT,Delta,a,tau)
[~,p] = size(Z);

alph = copulaparam('Frank',tau);  % %  the scalar parameter alpha.

beta = ini_beta;
phi = ini_phi;

lambda0 = log(sum(Delta));  % According to Su et al., 2016

lambda1 = 1e+7*sqrt(n);
lambda2 = 1e+7*sqrt(n);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1; err1=0; err2=0; err3=0; err4=0; 

tk1= 500; tk2= 500;

while k<=1000 && (err1==0 ||err2==0 || err3==0 || err4==0)
    %k
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

    DL_beta =  -(sum( (Delta-Delta.*H_T-LAM1.*S_T.*H_T).*Z,1) )'/n;  % % first partial derivative - beta
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
    uu = 0.6; % % taking value between 0 and 1
    gamm = 0.5; % % taking value between 0 and 0.5
    sigm = 0.35; % % taking value between 0 and 1

    t = 1; % % marking the loop

    while (t>0)

        H_T = sum( exp(Z*beta).*(seta.*Psi),2 );  % % Cumulative Risk - Failure
        H_C = sum( exp(Z*phi).*(gama.*Psi),2 );   % % Cumulative Risk - Censoring

        S_T = exp(-H_T);  % % Survival Probility - Failure
        S_C = exp(-H_C);  % % Survival Probility - Censoring

        S1 = exp(alph*S_T).*exp(alph*S_C)-exp(alph*S_T)-exp(alph*S_C)+exp(alph);  % % Copula: K(a,b,alph)
        S2 = exp(alph)*exp(alph*S_C)+exp(alph*S_C)-exp(2*alph*S_C)-exp(alph);
        S3 = exp(alph)*exp(alph*S_T)+exp(alph*S_T)-exp(2*alph*S_T)-exp(alph);

        DL_seta =  sum( (Delta.*psix)./(sum(seta.*psix,2))- ...
            (exp(Z*beta).*(Delta+LAM1.*S_T)).*Psi,1)/n - lambda1*seta*RT;   % % first partial derivative - seta

        DL_gama =  sum( ((1-Delta).*psix)./(sum(seta.*psix,2))- ...
            (exp(Z*phi).*(1-Delta+LAM2.*S_C)).*Psi,1)/n - lambda2*gama*RT;  % % first partial derivative - gama

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
            t = 0; % % break the loop
            uu_armijo = uu;
        else
            uu = uu*sigm; % %  reducing uu, enter the next loop
        end
    end

    seta1 = seta + uu_armijo*DL_seta*diag(seta./Ps1);
    gama1 = gama + uu_armijo*DL_gama*diag(gama./Ps2);
    err3 = norm(seta1-seta,2)^2 <= 1e-5*norm(seta,2)^2;
    err4 = norm(gama1-gama,2)^2 <= 1e-5*norm(gama,2)^2;
    gama = gama1;
    seta = seta1;

    k = k+1;
end

%%
beta2 = beta1.*(abs(beta1)>=8*1e-4);
phi2 = phi1.*(abs(phi1)>=8*1e-4);

opt_beta = beta2;
opt_phi = phi2;
opt_seta = seta;
opt_gama = gama;

end
