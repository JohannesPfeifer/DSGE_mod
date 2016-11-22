function [ys,check] = Jermann_Quadrini_2012_RBC_steadystate(ys,exo)
% [ys,check] = Jermann_Quadrini_2012_RBC_steadystate(ys,exo)
% computes the steady state for the Jermann_Quadrini_2012_RBC.mod and uses a numerical
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
xi_0=0.16338; %starting value
[xi,fval,exitflag]=fsolve(@(x)get_xi(x,tau,betta,delta,theta,BY_ratio),xi_0,options); %get xi that satisfies calibration targets
if exitflag<1
    check=1;
end
z=1;
n=0.3;
R=(1-tau)/betta+tau;
mu=(1-R*betta)/(xi*(R*(1-tau)/(R-tau)));
k=((((1-xi*mu)/betta)-(1-delta))/((1-mu)*theta*z*n^(1-theta)))^(1/(theta-1));
w=((1-theta)*k^theta*n^(-theta)*z)*(1-mu);
b=(z*k^theta*n^(1-theta)/xi-k)*(tau-R)/(1-tau);
d=(1-delta)*k+z*k^(theta)*n^(1-theta)-w*n-b+b/R-k;
c=w*n+b-b/R+d;
alppha=w/c^siggma*(1-n);

y=z*k^theta*n^(1-theta);
invest=delta*k;
v=d/(1-betta);
r=(R-tau)/(1-tau)-1;
yhat=0;
chat=0;
ihat=0;
nhat=0;
muhat=0;
byhat=0;
dyhat=d/y;
vyhat=0;


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


function resid=get_xi(xi,tau,betta,delta,theta,BY_ratio)
% compute deviation of debt-to-GDP ratio from target for given xi
z=1;
n=0.3;
R=(1-tau)/betta+tau;
r=(R-tau)/(1-tau)-1;
mu=(1-R*betta)/(xi*(R*(1-tau)/(R-tau)));
k=((((1-xi*mu)/betta)-(1-delta))/((1-mu)*theta*z*n^(1-theta)))^(1/(theta-1));
b=(z*k^theta*n^(1-theta)/xi-k)*(tau-R)/(1-tau);
y=z*k^theta*n^(1-theta);

resid=(BY_ratio-(b/(1+r)/y));

end