function resid=distance(shock_values,shock_name,target_value,target_name,M_,options_,oo_)
%function resid=distance(shock_values,shock_name,target_value,target_name,M_,options_,oo_)
% computes distance between target path and actual simulated path
% Inputs:
% - shock_values    [double]    t by 1 vector of shock values 
% - shock_name      [char]      name of the shock
% - target_value    [double]    t by 1 vector of target values 
% - target_name     [char]      name of the target variable
% - M_              [struct]    Dynare Model structure
% - options_        [struct]    Dynare options structure
% - oo_             [struct]    Dynare results structure
% 
% Outputs:
% - resid           [double]    distance vector
%

% Copyright (C) 2021 Johannes Pfeifer
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

oo_.exo_simul(M_.maximum_lag+1:M_.maximum_lag+length(shock_values),strmatch(shock_name,M_.exo_names,'exact'))=shock_values;
[oo_.endo_simul, success]= perfect_foresight_solver_core(oo_.endo_simul,oo_.exo_simul,oo_.steady_state,oo_.exo_steady_state,M_,options_);
if ~success
    resid=repmat(options_.huge_number,size(target_value));
else
    resid=(target_value-oo_.endo_simul(strmatch(target_name,M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+length(target_value)))';
end