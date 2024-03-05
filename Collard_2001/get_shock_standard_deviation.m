% This file shows how to find the shock standard deviation generating a
% specified resonse to a variable in a given period. It relies on internal 
% Dynare routines to make sure that specified correlations are preserved
% and the covariance matrix stays positive definite

% Copyright (C) 2023 Johannes Pfeifer
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

dynare Collard_2001_example1; %run Dynare once on the model to set up relevant options.
 
%% Options to be set by user
target_variable='y';  
shock_name='e';
target_value=1; % Desired deviation from steady-state in target period
target_period=2; % Period in which target should be achieved (impact period: 1)
starting_value=0.1; % Starting value of standard deviation for solver

%% Implementation

%set relevant estim_params_-structure
estim_params_.nvx=1;
estim_params_.ncx=0;
estim_params_.nvn=0;
estim_params_.ncn=0;
estim_params_.np=0;
estim_params_.var_exo=NaN(1,9);
estim_params_.var_exo(1,1)=strmatch(shock_name,M_.exo_names,'exact');
bounds.lb=0;
bounds.ub=Inf;
estim_params=check_for_calibrated_covariances(estim_params_,M_);

%run solver
options_.noprint=1;
options_.nomoments=1;

[x,exitflag] = csolve('minimize_distance',starting_value,[],1e-6,100,M_,options_,oo_,estim_params_,bounds,target_value,target_period,target_variable,shock_name);

%display info
if exitflag==0
    M_ = set_all_parameters(x,estim_params_,M_);
    [info, oo_] = stoch_simul(M_, options_, oo_, {target_variable});
    fprintf('\nShock standard deviation of %1.3f will yield a response of %1.3f to variable %s in period %d.\n',x,oo_.irfs.([target_variable '_' shock_name])(target_period), target_variable, target_period)
else
    error('No solution could be found')
end