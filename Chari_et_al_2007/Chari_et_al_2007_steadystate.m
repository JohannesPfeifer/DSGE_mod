function [ys,check] = Chari_et_al_2007_steadystate(ys,exo)
% function [ys,check] = Chari_et_al_2007_steadystate(ys,exo)
% computes the steady state for the Chari_et_al_2007.mod; when doing so,
% the estimated coefficients are translated into the VAR representation for
% the exogenous processes
% 
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
%construct P0-matrix
P0=[P0_z_bar;P0_tau_l_bar;P0_tau_x_bar;P0_g_bar]; 
%construct P-matrix
P=[rho_zz rho_zl rho_zx rho_zg
    rho_lz rho_ll rho_lx rho_lg
    rho_xz rho_xl rho_xx rho_xg
    rho_gz rho_gl rho_gx rho_gg];
%construct \bar s vector of means
s_bar=(eye(4,4)-P)\P0;
%read out means
z_bar=s_bar(1);
tau_l_bar=s_bar(2);
tau_x_bar=s_bar(3);
g_bar=s_bar(4);

tau_l=tau_l_bar;
z=z_bar;
tau_x=tau_x_bar;
g=g_bar;
beta_hat=betta*(1+gamma_z)^(-siggma);
cap_labor_ratio=(((1+tau_x)*(1-beta_hat*(1-delta)))/(beta_hat*theta*exp(z)^(1-theta)))^(1/(theta-1));
xi_1=cap_labor_ratio^(theta-1)*exp(z)^(1-theta)-(1+gamma_z)*(1+gamma_n)+(1-delta);
xi_2=(1-tau_l)*(1-theta)*cap_labor_ratio^theta*exp(z)^(1-theta)/psii;
xi_3=xi_2/cap_labor_ratio;
k=(xi_2+exp(g))/(xi_1+xi_3);
c=xi_1*k-exp(g);
l=1/cap_labor_ratio*k;
x=(1+gamma_n)*(1+gamma_z)*k-(1-delta)*k;
w=(1-theta)*k^theta*l^(-theta)*exp(z)^(1-theta);
y=k^theta*(exp(z)*l)^(1-theta);

log_labor_wedge = log(1-tau_l);
log_investment_wedge = log(1/(1+tau_x));
log_efficiency_wedge=log(exp(z)^(1-theta));

k=log(k);
w=log(w);
c=log(c);
y=log(y);
x=log(x);
l=log(l);

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
