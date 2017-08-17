function W=return_welfare(xopt)
% function W=return_welfare(xopt)
% Sets parameter and returns objective function

% Copyright (C) 2017 Johannes Pfeifer
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

global oo_ M_ options_ % get Dynare structures; used to pass them on to resol.m

%% set parameter for use in Dynare
set_param_value('tau_n',xopt(1));

info=stoch_simul('W'); %run computation of mean welfare

if info(1) %solution was not successful
    W=10e6+sum([xopt(1)].^2); %return with penalty 
else
    W=-oo_.mean(1); %return minus welfare as we are using minimizer
end
