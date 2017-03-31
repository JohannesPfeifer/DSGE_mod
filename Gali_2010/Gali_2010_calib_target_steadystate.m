function [ys,check] = Gali_2010_calib_target_steadystate(ys,exo)
% function [ys,check] = Gali_2010_calib_target_steadystate(ys,exo)
% uses a numerical solver to get the steady state for the nonlinear version of Gali_2010_calib_target_steadystate.mod 
% in order to set the required parameters matching the calibration
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
% Copyright (C) 2016 Johannes Pfeifer
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% For a copy of the GNU General Public License,
% see <http://www.gnu.org/licenses/>.

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

% definition labor force, p. 516
F=N+U;

% steady state separation rate, p.516
delta=x/(1-x)*U/N;        


%% solve (48) to (52) for psi, chi, C, L numerically
if matlab_ver_less_than('9.1')
    options=optimset('TolFun',1e-10,'TolX',1e-10);
else
    options=optimoptions('fsolve','FunctionTolerance',1e-15,'OptimalityTolerance',1e-15);
end
[outvalue,fval,exitflag]=fsolve(@(par_vec)solve_calibration(par_vec,betta,delta,hiring_cost_share,x,alfa,N,xi,U,varphi),[0.041;15.5;0.7;0.05],options);
if exitflag <1
    %indicate the SS computation was not sucessful; this would also be detected by Dynare
    %setting the indicator here shows how to use this functionality to
    %filter out parameter draws
    check=1; %set failure indicator
    return; %return without updating steady states
end

psi=outvalue(1);
chi=outvalue(2);
G=outvalue(3);
C=outvalue(4);

%p. 499 bottom in steady state 
Gamma=G/x^gammma;

MRPN=(1-alfa)*N^(-alfa);
%equation (51), p. 515
L=N+psi*U;
MRS   = chi*C*L^varphi;
W_div_P    = xi*MRS + (1-xi)*MRPN; % 
Upsilon=xi*MRS/W_div_P; %p. 510 below (35)
B=(1-betta*(1-delta))*G;
Phi=B/(W_div_P+B);% p. 502
Y=N^(1-alfa);
%Coefficient optimal participation condition, bottom p. 541 
Xi=(xi*W_div_P/(G*(1-xi)))*(theta_w/((1-theta_w)*(1-betta*theta_w*(1-delta))));
Theta=delta*N*G/Y;
fprintf('The ratio of average hiring costs to quarterly wage is %4.3f\n',G/W_div_P)
%% The conditions that need to be satisfied
residuals=NaN(6,1);
residuals(1)=N^(1-alfa)-(C+delta*N*Gamma*x^gammma); %(48)
residuals(2)=(1-betta*(1-delta))*Gamma*x^gammma-xi*((1-alfa)*N^(-alfa)-chi*C*L^varphi);%(49)
residuals(3)=(1-x)*xi*psi*chi*C*L^varphi-(1-xi)*Gamma*x^(1+gammma); %(50)
residuals(4)=x*U-(1-x)*delta*N; %(51)
residuals(5)=L-(N+psi*U); %(52)
residuals(6)=G/W_div_P-hiring_cost_share;

if max(abs(residuals)>1e-6)
    check=1; %set failure indicator
    return; %return without updating steady states
end


%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = 0;']);
end
end