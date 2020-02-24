function outvalue=get_consumption_equivalent_unconditional_welfare(par_value_lambda)
% function outvalue=get_consumption_equivalent_unconditional_welfare(par_value)
% computes the unconditional welfare difference between the given model equilibrium and the 
% alternative with a consumption-equivalent of lambda lost compared to the
% natural/flex-price allocation; it relies on the unconditional mean of
% lifetime utility defined in the model and computed by Dynare

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

global oo_ options_ M_

if par_value_lambda>1 %do not allow negative consumption
    outvalue=1e5+par_value_lambda^2;
end

set_param_value('lambda_utility',par_value_lambda) %set consumption equivalent lambda
if ~options_.ramsey_policy
    var_list_ = {'Welfare_gap'};
    [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %get decision rules and moments
    if info(1) %filter out error code
        outvalue=1e5+par_value_lambda^2;
        return;
    end
    outvalue=oo_.mean(strmatch('Welfare_gap',var_list_,'exact')); %extract Welfare gap measure
else
    %workaround for Ramsey where recursive welfare cannot be defined in the
    %model; uses the fact that E(W)=E(U)+betta*E(W) defines the fixed point
    %for welfare
    var_list_ = char('Utility','Recursive_natural_welfare_equivalent');
    info = stoch_simul(var_list_);
    if info(1)
        outvalue=1e5+par_value_lambda^2;
        return;
    end
    betta=M_.params(strmatch('betta',M_.param_names,'exact'));
    outvalue=1/(1-betta)*oo_.mean(strmatch('Utility',var_list_ ,'exact'))-oo_.mean(strmatch('Recursive_natural_welfare_equivalent',var_list_,'exact'));
end