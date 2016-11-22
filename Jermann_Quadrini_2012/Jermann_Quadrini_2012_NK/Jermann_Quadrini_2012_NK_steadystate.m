function [ys,check] = Jermann_Quadrini_2012_NK_steadystate(ys,exo)
% function [ys,check] = Jermann_Quadrini_2012_NK_steadystate(ys,exo)
% computes the steady state for the Jermann_Quadrini_2012_NK.mod and uses a numerical
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
%
% Copyright (C) 2014-2016 Johannes Pfeifer
% 
%  This is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  It is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  For a copy of the GNU General Public License,
%  see <http://www.gnu.org/licenses/>.

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

G=G_bar;
z=0;
zeta=1;
gamma=1;
eta=eta_bar;
var_sigma=0;
upsilon=upsilon_bar;

P=1;
Q=1;
chi=0;
u=1;

xi_0=0.2;
[xi,fval,exitflag]=fsolve(@(x)get_xi(x,tau,betta,delta,theta,eta,BY_ratio),xi_0,options); %get xi that satisfies calibration targets
xi_bar=xi;

if exitflag~=1
    check=1;
    return
end
n=0.3;
R=(1-tau)/betta+tau;
r=(R-tau)/(1-tau)-1;
mu=(1-R*betta)/(xi*R/(1+r));

k=((((1-xi*mu)/betta)-(1-delta))/((1-mu)*theta/eta*n^(1-theta)))^(1/(theta-1));
W=1/eta*(1-theta)*k^theta*n^(-theta)*(1-mu);
w_opt=W;
b=(k^theta*n^(1-theta)/xi-k)*(tau-R)/(1-tau);
y=k^theta*n^(1-theta);
G=GY_ratio*y;
G_bar=G;
T=G+b*(1/R-1/(1+r));
d=(1-delta)*k+k^(theta)*n^(1-theta)-W*n-b+b/R-k;
c=W*n+b-b/(1+r)+d-T;
alppha=W/(upsilon_bar*(c*(1-h))^siggma*n^(1/epsilon)); %using that w=upsilon*MRS=upsilon*n^(1/epsilon)/((c*(1-h))^-siggma)

invest=delta*k;
V=d/(1-betta);
yhat=0;
chat=0;
ihat=0;
nhat=0;
muhat=0;
byhat=0;
dyhat=d/y;
vyhat=0;
wopthat=0;
upsilonhat=0;
What=0;

y_obs =0;
c_obs =0;
invest_obs  =0;
pi_obs  =0;
r_obs  =0;
n_obs  =0;
W_obs  =0;
debt_repurchase_obs =0;

%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end

end


function resid=get_xi(xi,tau,betta,delta,theta,eta,BY_ratio)
% compute deviation of debt-to-GDP ratio from target for given xi

n=0.3;
R=(1-tau)/betta+tau;
r=(R-tau)/(1-tau)-1;
mu=(1-R*betta)/(xi*R/(1+r));

k=((((1-xi*mu)/betta)-(1-delta))/((1-mu)*theta/eta*n^(1-theta)))^(1/(theta-1));
b=(k^theta*n^(1-theta)/xi-k)*(tau-R)/(1-tau);
y=k^theta*n^(1-theta);

resid=(BY_ratio-(b/(1+r)/y));

end

