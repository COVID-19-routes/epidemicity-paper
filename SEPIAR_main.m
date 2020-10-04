clearvars
close all
clc

%% setup
% population distribution and mobility
load data_provinces pop mob_pop mob_mat
N=pop; % resident population in each province
p=mob_pop; % fraction of mobile people in each province
q=mob_mat; % province-to-province mobility matrix
% data_provinces.mat also contains:
% id: ISTAT code of each province
% x,y: longitude and latitude of provincial centroids
clear pop mob_pop mob_mat
n=numel(N); % number of provinces

% parameter values
% Bertuzzo et al. (2020), Nat. Comm. 11:4264
% https://doi.org/10.1038/s41467-020-18050-2
mu=1/75/365; % baseline mortality rate (1/days)
deltaE=1/4.6; % latency rate (1/days)
deltaP=1/2; % post-latency rate (1/days)
sigma=0.25; % fraction of symptomatics
eta=1/5; % hospitalization rate of symptomatics (1/days)
gammaI=1/14; % recovery rate of symptomatics (1/days)
alpha=1/25; % mortality rate of symptomatics (1/days)
gammaA=2*gammaI; % recovery rate of asymptomatics (1/days)
r_hi=0.5; % fraction of mobile contacts for S/E/P/A/R
r_lo=0; % fraction of mobile contacts for I
rS=r_hi; rE=r_hi; rP=r_hi; rI=r_lo; rA=r_hi; rR=r_hi;
betaP=1.26; % transmission rate from post-latent (1/days)
betaI=0.022*betaP; % transmission rate from symptomatics (1/days)
betaA=betaI; % transmission rate from asymptomatics (1/days)

% some useful matrices
zv=zeros(n,1); % all-zero column vector of size n
Zm=sparse(n,n); % all-zero n-by-n matrix
uv=ones(n,1); % all-one column vector of size n
Um=ones(n,n); % all-one n-by-n matrix
u1n=1:n; % 1-to-n row vector
U=sparse(u1n,u1n,uv,n,n); % identity matrix of size n

% output matrix for the epidemicity index
% (all and only the state variables of the infection subsystem have
% weight = 1 in the output transformation)
W=diag(ones(10*n,1)); W([u1n n*5+1:end],:)=[];
Wpi=sparse(pinv(W)); % pseudoinverse of W
W=sparse(W);

% simulation presettings
ndays=60; % number of days in the simulation
tspan=0:ndays;
prov=15; % province with the first transmission focus (e.g. Milano)
E0=zv; E0(prov)=10; % number of exposed in the first transmission focus (e.g. 10)
x0=[N-E0; E0; repmat(zv,8,1); E0];

%% without controls
epsilon=zv;
xi=Zm;
chi=zv;
chiE=chi; chiP=chi; chiI=chi; chiA=chi;

% contact matrices
p_diag=sparse(u1n,u1n,p,n,n); clear p
MS=rS*p_diag*(q.*(Um-xi))+(1-rS)*p_diag+U-p_diag;
ME=rE*p_diag*(q.*(Um-xi))+(1-rE)*p_diag+U-p_diag;
MP=rP*p_diag*(q.*(Um-xi))+(1-rP)*p_diag+U-p_diag;
MI=rI*p_diag*(q.*(Um-xi))+(1-rI)*p_diag+U-p_diag;
MA=rA*p_diag*(q.*(Um-xi))+(1-rA)*p_diag+U-p_diag;
MR=rR*p_diag*(q.*(Um-xi))+(1-rR)*p_diag+U-p_diag;

% basic reproduction number and epidemicity index
% (effective reproduction number/epidemicity index at time 0 without controls)
[R0,e0]=SEPIAR_eigen_t(...
    betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
    mu,deltaE,deltaP,sigma,gammaI,eta,alpha,gammaA,...
    N,zv,zv,zv,zv,zv,MS,ME,MP,MI,MA,MR,n,W,Wpi);

disp('***without controls***')
disp(['basic reproduction number = ',num2str(R0)])
disp(['basic epidemicity index = ',num2str(e0),' (1/days)'])

% simulation
ode=@(t,x) SEPIAR_ode(t,x,...
    betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
    mu,deltaE,deltaP,sigma,gammaI,eta,alpha,gammaA,...
    N,MS,ME,MP,MI,MA,MR,n);
[~,x]=ode45(ode,tspan,x0);

St=x(:,u1n); % susceptible people, (ndays+1)xn matrix
Et=x(:,n+1:2*n); % exposed people
Pt=x(:,2*n+1:3*n); % post-latent people
It=x(:,3*n+1:4*n); % symptomatic people
At=x(:,4*n+1:5*n); % asymptomatic people
Rt=x(:,9*n+1:10*n); % recovered people
Ct=x(:,10*n+1:end); % cumulated cases

% total cumulated cases and total prevalence at the end of the simulation
totcases=sum(Ct(end,:));
totprev=sum(Et(end,:)+Pt(end,:)+It(end,:)+At(end,:))/...
    sum(St(end,:)+Et(end,:)+Pt(end,:)+It(end,:)+At(end,:)+Rt(end,:));

disp(['total cumulated cases = ',num2str(totcases/1e3),' (thousands)'])
disp(['total prevalence = ',num2str(totprev*1e5),' (per 100 thousand population)'])

% total daily cases and cumulated cases by province
tdc_nc=[sum(Ct(1,:)) diff(sum(Ct,2))'];
ccp_nc=Ct(end,:);

%% with controls
epsilon=uv*.3;
xi=Um*.5;
chi=uv*.1; % (1/days)
chiE=chi; chiP=chi; chiI=chi; chiA=chi;

% contact matrices
MS=rS*p_diag*(q.*(Um-xi))+(1-rS)*p_diag+U-p_diag;
ME=rE*p_diag*(q.*(Um-xi))+(1-rE)*p_diag+U-p_diag;
MP=rP*p_diag*(q.*(Um-xi))+(1-rP)*p_diag+U-p_diag;
MI=rI*p_diag*(q.*(Um-xi))+(1-rI)*p_diag+U-p_diag;
MA=rA*p_diag*(q.*(Um-xi))+(1-rA)*p_diag+U-p_diag;
MR=rR*p_diag*(q.*(Um-xi))+(1-rR)*p_diag+U-p_diag;

% control reproduction number and epidemicity index
% (effective reproduction number/epidemicity index at time 0 with controls)
[Rc,ec]=SEPIAR_eigen_t(...
    betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
    mu,deltaE,deltaP,sigma,gammaI,eta,alpha,gammaA,...
    N,zv,zv,zv,zv,zv,MS,ME,MP,MI,MA,MR,n,W,Wpi);

disp('***with controls***')
disp(['control reproduction number = ',num2str(Rc)])
disp(['control epidemicity index = ',num2str(ec),' (1/days)'])

% simulation
ode=@(t,x) SEPIAR_ode(t,x,...
    betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
    mu,deltaE,deltaP,sigma,gammaI,eta,alpha,gammaA,...
    N,MS,ME,MP,MI,MA,MR,n);
[~,x]=ode45(ode,tspan,x0);

St=x(:,u1n); % susceptible people, (ndays+1)xn matrix
Et=x(:,n+1:2*n); % exposed people
Pt=x(:,2*n+1:3*n); % post-latent people
It=x(:,3*n+1:4*n); % symptomatic people
At=x(:,4*n+1:5*n); % asymptomatic people
Rt=x(:,9*n+1:10*n); % recovered people
Ct=x(:,10*n+1:end); % cumulated cases

% total cumulated cases and total prevalence at the end of the simulation
totcases=sum(Ct(end,:));
totprev=sum(Et(end,:)+Pt(end,:)+It(end,:)+At(end,:))/...
    sum(St(end,:)+Et(end,:)+Pt(end,:)+It(end,:)+At(end,:)+Rt(end,:));

disp(['total cumulated cases = ',num2str(totcases/1e3),' (thousands)'])
disp(['total prevalence = ',num2str(totprev*1e5),' (per 100 thousand population)'])

% total daily cases and cumulated cases by province
tdc_wc=[sum(Ct(1,:)) diff(sum(Ct,2))'];
ccp_wc=Ct(end,:);

%% comparison
figure
subplot(211)
plot(tspan,tdc_nc/1e3,tspan,tdc_wc/1e3,'LineWidth',1)
xlabel('time t (days)')
ylabel('total daily cases (thousands)')
legend('without controls','with controls','Location','northwest')
subplot(212)
bar(u1n,[ccp_nc; ccp_wc]',1)
set(gca,'YScale','log','YLim',[1 1e6])
xlabel('provinces by increasing id')
ylabel('cumulated cases')
