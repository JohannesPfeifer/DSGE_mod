function output_cell = posterior_variance_decomposition(xparam1,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info),
% function output_cell = posterior_variance_decomposition(xparam1,M_,options_,oo_,estim_params_,bayestopt_,dataset_,dataset_info),
% Function called by execute_prior_posterior_function that generates
% variance decomposition for given set of parameters
% For input/output arguments, see the Dynare manual on posterior_function

% Copyright (C) 2016 Johannes Pfeifer
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

var_list=char('y_obs','c_obs','invest_obs','pi_obs','r_obs','n_obs','W_obs','debt_repurchase_obs');
M_ = set_all_parameters(xparam1,estim_params_,M_);
options_.noprint=1;
[dr,info,M_local,options_local,oo_local] = resol(0,M_,options_,oo_);
if ~info
    oo_local=disp_th_moments(dr,var_list,M_local,options_local,oo_local);
    output_cell={oo_local.variance_decomposition,xparam1};
else
    output_cell={[],[]};
end