function outvalue=get_consumption_equivalent_conditional_welfare(par_value_lambda)
% function outvalue=get_consumption_equivalent_conditional_welfare(par_value_lambda)
% computes the conditional welfare difference between the given model equilibrium and the 
% alternative with a consumption-equivalent of lambda lost compared to the
% natural/flex-price allocation; it uses Dynare's simult_-function to compute the lifetime utility 
% defined in the model for a given state vector (here the steady state)

% Copyright (C) 2018 Johannes Pfeifer and Benjamin Born
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

global oo_ M_ options_

if par_value_lambda>1 %do not allow negative consumption
    outvalue=1e5+par_value_lambda^2;
end

set_param_value('lambda_utility',par_value_lambda)  %set consumption equivalent lambda
if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end
[oo_.dr,info,M_,options_,oo_] = resol(0,M_,options_,oo_); %get decision rules
if info(1) %filter out error codes
    outvalue=1e5+par_value_lambda^2;
    return;
end

%% simulate conditional welfare
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag); %get steady state as initial condition
shock_matrix = zeros(1,M_.exo_nbr); %create shock matrix with number of time periods in rows
y_sim = simult_(M_,options_,initial_condition_states,oo_.dr,shock_matrix,options_.order); %simulate one period to get value

if ~options_.ramsey_policy
    outvalue=y_sim(strmatch('Welfare_gap',M_.endo_names,'exact'),2)*100; %read out gap
else
    Recursive_natural_welfare_equivalent=y_sim(strmatch('Recursive_natural_welfare_equivalent',M_.endo_names,'exact'),2);
    oo_.steady_state=oo_.dr.ys;
    oo_.planner_objective_value = evaluate_planner_objective(M_,options_,oo_);
    outvalue=(oo_.planner_objective_value(1)-Recursive_natural_welfare_equivalent)*100;
end