function distance=minimize_distance(sigma,M_,options_,oo_,estim_params_,bounds,...
                                target_value,target_period,target_variable,shock_name)
%function outvalue=minimize_distance(sigma,M_,options_,oo_,estim_params_,bounds,...
%                                target_value,target_period,target_variable,shock_name)
%
% Input:
% - sigma           [double]    value of standard deviation
% - M_              [struct]    Dynare model structure
% - options_        [struct]    Dynare options structure
% - oo_             [struct]    Dynare results structure 
% - estim_params_   [struct]    Dynare estimated parameter information structure
% - bounds          [struct]    Dynare structure specifying bounds on the
%                               estimated parameter (subfields lb and ub)
% - target_value    [double]    target value for the IRF
% - target_period   [integer]   period in which target should be achieved
% - target_variable [char]      name of target variable
% - shock_name      [char]      name of the shock for which to set standard
%                               deviation
%
% Output:
% -distance     [double]    distance from target 

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

M_ = set_all_parameters(sigma,estim_params_,M_); %set parameter vector

%check parameter bounds and definiteness
try
    [fval,info]=check_bounds_and_definiteness_estimation(sigma, M_, estim_params_, bounds);
catch % Dynare versions < 5.0
    fval        = [];
    exit_flag   = 1;
    info        = zeros(4,1);
    % Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
    if sigma<bounds.lb
        info(1) = 41;
        info(4) = sum((bounds.lb-sigma).^2);
    end
    % Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
    if sigma>bounds.ub        
        info(1) = 42;
        info(4) = sum((sigma-bounds.ub).^2);
    end
end
if info(1)
    distance=options_.huge_number;    
    return
end
          
%get IRF and compute distance
[info, oo_] = stoch_simul(M_, options_, oo_, {target_variable});
if info~=0
    distance=options_.huge_number;
else
    distance=oo_.irfs.([target_variable '_' shock_name])(target_period)-target_value;
end