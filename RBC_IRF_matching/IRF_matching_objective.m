function [fval, IRF_model]=IRF_matching_objective(xopt,IRF_target,weighting_matrix)
% function  [fval, IRF_model]=IRF_matching_objective(xopt,IRF_target,weighting_matrix)
% Computes the quadratic deviation of the simulated IRFs from the empirical
% ones; does so by calling stoch_simul.m from Dynare to compute IRFs
%
% IRFs
% Inputs:
%   xopt                [npar by 1]             vector of parameters in current optimization step
%   IRF_target          [nperiods by nvars]     matrix of target IRFs
%   weighting_matrix    [nperiods*nvars by nperiods*nvars]     matrix of weights for IRF matching 
%                                                               (often a matrix with the inverse of the 
%                                                               variances of the IRFS on the diagonal)
% Outputs:
%   fval            [scalar]                value of the objective function (quadratic distance)
%   IRF_model       [nperiods by nvars]     model IRFs
% 
% Notes:
%   - when calling stoch_simul.m within a loop/optimizer, it is essential
%       that options_.noprint has been set. This will prevent Dynare from
%       throwing errors and allows using the error code returned by
%       stoch_simul to add a penalty
%   - The empirical IRFs and model IRFs use an impulse size of 1 percent. Thus, there is no uncertainty about the
%       initial impact. The IRF matching therefore only targets the G-response starting in the second period.

% Copyright (C) 2017 Johannes Pfeifer
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
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.

global oo_ M_ options_ %required input structures for call to stoch_simul

%% set parameter for use in Dynare
set_param_value('root_g_1',xopt(1)); % first autoregressive root
set_param_value('root_g_2',xopt(2)); % second autoregressive root

if any(xopt<=-1) || any(xopt(1)>=1) %make sure roots are between 0 and 1
    fval=10e6+sum([xopt].^2); %penalty function
    return
end

info=stoch_simul([]); %run stoch_simul to generate IRFs with the options specified in the mod-file

if info %solution was not successful
    fval=10e6+sum([xopt].^2); %return with penalty 
else
    IRF_model=[oo_.irfs.ghat_eps_g' oo_.irfs.log_y_eps_g']; %select desired IRFs
    % compute objective function (omitting the impact response of G as that is targeted with the shock size)
    fval=(IRF_model(2:end)-IRF_target(2:end))*weighting_matrix*(IRF_model(2:end)-IRF_target(2:end))';  %start only in the second period for g 
end