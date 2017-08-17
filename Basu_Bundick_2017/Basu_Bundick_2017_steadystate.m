function [ys,check] =BB2016_steadystate(ys,exe)
% function [ys,check] = NK_baseline_steadystate(ys,exo)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values fpr the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impos restriction on parameters)

global M_ 

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% Enter model equations here
theta_v =(1-siggma)/(1-1/psii);
V=1;
Y=1;
PHI=1;
M=betta;
E_t_V_tp1_1_minus_sigma=1;
Xi = (theta_mu-1)/theta_mu;
mu=1/Xi;
Pi=Pi_bar;
R=Pi_bar/betta;
q=1;
R_K=1/betta-(1-delta_0);
S=1;
a=a_bar;
Z=1;
sigma_a=sigma_a_bar;
R_R=1/betta;
delta_1 = 1/betta-1+delta_0;
u = 1;

K=alppha/R_K;
options=optimset('TolFun',10e-12,'TolX',10e-12);
%N=fsolve(@(N) (1-alppha)*(1-N-siggma*N)/(siggma*((2*N-1)/(1-N))*N)+delta*K-1,0.4)
N=fsolve(@(N) (1-alppha)/N*(1-N)*(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N)))/(1-(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N))))+delta_0*K-1,0.02,options)
eta=(theta_v/(1-siggma)*(1-Frisch_target*(1-(1-siggma)/theta_v)*N/(1-N)));
Frisch=(1-eta*(1-siggma)/theta_v)/(1-(1-siggma)/theta_v)*(1-N)/N


% N=fsolve(@(N) (1-alppha)/N*(1-N)*(1/(1-siggma)-siggma/(1-siggma)*N/(1-N))/(1-(1/(1-siggma)-siggma/(1-siggma)*N/(1-N)))+delta*K-1,0.4,options);
% eta=1/(1-siggma)-siggma/(1-siggma)*N/(1-N);
% Frisch=(1-eta*(1-siggma))/siggma*(1-N)/N;
if abs(Frisch-Frisch_target)>10e-5
   error('Wrong calibration. Frisch elasticity is not at target') 
end


PF_normalization=(1/(Xi*N^(1-alppha)*K^alppha));
Y_gross=Xi*(PF_normalization*N^(1-alppha)*K^alppha);
W=(1-alppha)*Y_gross/N;
Psi=(1-Xi)*PF_normalization*N^(1-alppha)*K^alppha;
D=Y-N*W-delta_0*K;
B=nu*K;
D_E=D-(1-1/R_R)*B;
P_E=betta/(1-betta)*D_E;
I=delta_0*K;
C=Y-I;

V_normalization=(1-betta)/((C^eta*(1-N)^(1-eta))^((1-siggma)/theta_v));
R_E=(D_E + P_E)/P_E;
E_R_E=R_E;
cond_var_R_E=0;
E_R_E_squared=R_E^2;

E_M_tp1=betta;
E_R_E_risk_neutral=R_E/E_M_tp1;
E_M_tp1_R_E=1;
E_M_tp1_R_E_squared=betta*R_E^2;
cond_var_R_E_risk_neutral=E_M_tp1_R_E_squared/E_M_tp1-(E_M_tp1_R_E/E_M_tp1)^2;
E_RN_R_E=R_E;
E_RN_R_E_squared=R_E^2;
conditional_variance_R_E_risk_neutral_2=E_RN_R_E_squared-E_RN_R_E^2;

log_Y=log(Y);
log_C=log(C);
log_I=log(I);
log_mu=log(mu);
log_N=log(N);
log_sigma_a=log(sigma_a);
log_W =log(W);
pi_annualized=1+(4*(Pi-1));
R_annualized=1+(4*(R-1));
R_K_annualized=4*R_K;
R_R_annualized=1+(4*(R_R-1));
log_pi_annualized=log(pi_annualized);
log_R_annualized=log(R_annualized);
log_R_K_annualized=log(R_K_annualized);
log_R_R_annualized=log(R_R_annualized);

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
