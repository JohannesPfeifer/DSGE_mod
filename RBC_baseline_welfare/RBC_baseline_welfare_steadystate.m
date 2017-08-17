function [ys,check] = RBC_baseline_welfare_steadystate(ys,exo)
% function [ys,check] = RBC_baseline_welfare_steadystate(ys,exo)
% computes the steady state for the RBC_baseline_welfare.mod and uses a numerical
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

options=optimset('display','off'); % set options for numerical solver

%Do Calibration
gammax=(1+n)*(1+x);
delta=i_y/k_y-x-n-n*x;
betta=(1+x)*(1+n)/(alppha/k_y+(1-delta));
[l,fval,exitflag]=fsolve(@(l)((1-tau_n)*(1-alppha)*(((1/betta*(1+n)*(1+x)-(1-delta))/alppha)^(alppha/(alppha-1)))...
    -psii*((((1/betta*(1+n)*(1+x)-(1-delta))/alppha)^(1/(alppha-1))*l)^(alppha)*l^(1-alppha)-(x+n+delta+n*x)*((1/betta*(1+n)*(1+x)-(1-delta))/alppha)^(1/(alppha-1))*l)^siggma/(1-l)),...
    0.25,options);
if exitflag <1
    %indicate the SS computation was not sucessful; this would also be detected by Dynare
    %setting the indicator here shows how to use this functionality to
    %filter out parameter draws
    check=1; %set failure indicator
    return; %return without updating steady states
end
k = ((1/betta*(1+n)*(1+x)-(1-delta))/alppha)^(1/(alppha-1))*l; 
invest = (x+n+delta+n*x)*k;
y=k^alppha*l^(1-alppha);
c = k^(alppha)*l^(1-alppha)-invest;
W=1/(1-betta)*(c^(1-siggma)/(1-siggma)+psii*log(1-l));
w = (1-alppha)*y/l;
r = 4*alppha*y/k;
log_y = log(y);
log_k = log(k);
log_c = log(c);
log_l = log(l);
log_w = log(w);
log_invest = log(invest);
z = 0; 

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
