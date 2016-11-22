function generate_var_decomp_posterior()
% function generate_var_decomp_posterior()
% Generates the prior/posterior variance decomposition graph
%
% Requirements: 
% - a previous estimation must have been conducted so that
%       Metropolis-files are located in the corresponding Dynare-folder
% - the kernel density estimator makes use of Dynare routines. To add them
%       to your path, you will have to run dynare_config.m if you have
%       not executed Dynare 4.5 in this Matlab session before calling this
%       function

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

global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info

load('Jermann_Quadrini_2012_NK_results.mat')
%reconstruct var_list from mod-file which was not saved
var_list=char('y_obs','c_obs','invest_obs','pi_obs','r_obs','n_obs','W_obs','debt_repurchase_obs');
xi_pos=strmatch('eps_xi',M_.exo_names,'exact');
y_pos=strmatch('y_obs',var_list,'exact');

% set up options required for variance decomposition
dataset_info=[];
options_.sampling_draws = 10000;
options_.use_qzdiv=1;
options_.function = 'posterior_variance_decomposition';

oo_ = execute_prior_posterior_function('posterior_variance_decomposition', M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, 'posterior');

var_decom_array_posterior=NaN(size(oo_.posterior_function_results{1,1},1),size(oo_.posterior_function_results{1,1},2),options_.sampling_draws);
for replic_iter=1:options_.sampling_draws
    var_decom_array_posterior(:,:,replic_iter)=oo_.posterior_function_results{replic_iter,1};
end

%plot results
figure('Name','Posterior')
plot(squeeze(var_decom_array_posterior(y_pos,xi_pos,:)))
mean(var_decom_array_posterior,3)
figure('Name','Posterior')
hist(squeeze(var_decom_array_posterior(y_pos,xi_pos,:)),40)

%display table
title='Posterior variance Decomposition (in percent)';
headers=M_.exo_names_tex;
headers = char(' ',headers);
labels = deblank(var_list);
lh = size(labels,2)+2;
dyn_latex_table(M_,title,'posterior_var_decomp_uncond',headers,labels,mean(var_decom_array_posterior,3),lh,8,2);

%% %%%%%%%%%%%%%%%%%% Do prior predictive plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oo_ = execute_prior_posterior_function('posterior_variance_decomposition', M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, 'prior');

var_decom_array_prior=NaN(size(oo_.posterior_function_results{1,1},1),size(oo_.posterior_function_results{1,1},2),options_.sampling_draws);
for replic_iter=1:options_.sampling_draws
    if ~isempty(oo_.prior_function_results{replic_iter,1})
        var_decom_array_prior(:,:,replic_iter)=oo_.prior_function_results{replic_iter,1};
    end
end

%plot results
figure('Name','Prior')
plot(squeeze(var_decom_array_prior(y_pos,xi_pos,:)))
nanmean(var_decom_array_prior,3)
figure('Name','Prior')
hist(squeeze(var_decom_array_prior(y_pos,xi_pos,:)),40)


%% get Kernel density estimator
%use Dynare defaults for Kerndel density estimator
number_of_grid_points = options_.estimation.moments_posterior_density.gridpoints;
bandwidth = options_.estimation.moments_posterior_density.bandwidth;
kernel_function = options_.estimation.moments_posterior_density.kernel_function;

%prior
draws=squeeze(var_decom_array_prior(y_pos,xi_pos,:));
draws=draws(~isnan(draws));
number_of_draws=length(draws);

optimal_bandwidth_prior = mh_optimal_bandwidth(draws,number_of_draws,bandwidth,kernel_function);
[density_prior(:,1),density_prior(:,2)] = kernel_density_estimate(draws,number_of_grid_points,...
    number_of_draws,optimal_bandwidth_prior,kernel_function);

%posterior
draws=squeeze(var_decom_array_posterior(y_pos,xi_pos,:));
draws=draws(~isnan(draws));
number_of_draws=length(draws);

optimal_bandwidth_posterior = mh_optimal_bandwidth(draws,number_of_draws,bandwidth,kernel_function);
[density_posterior(:,1),density_posterior(:,2)] = kernel_density_estimate(draws,number_of_grid_points,...
    number_of_draws,optimal_bandwidth_posterior,kernel_function);

figure('Name','Prior/Posterior Variance share')
plot(density_prior(:,1),density_prior(:,2),'k-',density_posterior(:,1),density_posterior(:,2),'r-.','LineWidth',1.5)
axis([0 50 0 0.17])
legend('Prior','Posterior')
xlabel('Variance Share')
ylabel('Density Estimate')

%get value at the mode
x_mode=load('Jermann_Quadrini_2012_NK_mode');
M_local = set_all_parameters(x_mode.xparam1,estim_params_,M_);
[dr,info,M_local,options_local,oo_local] = resol(0,M_local,options_,oo_);
oo_local=disp_th_moments(dr,var_list,M_local,options_local,oo_local);
hold on
plot([oo_local.variance_decomposition(y_pos,xi_pos) oo_local.variance_decomposition(y_pos,xi_pos)],[0 0.17],'g--')
plot([46.4 46.4],[0 0.17],'k:')
print('prior_posterior_variance_decomposition','-depsc2')