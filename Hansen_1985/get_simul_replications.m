function [sim_array]=get_simul_replications(M_,options_)
% function [sim_array]=get_simul_replications(M_,options_)
% reads the simulation replications into a three-dimensional array with
% endogenous variables along the first, simulation periods along the
% second, and replications along the third dimension
%
% INPUTS
%   M_:  Dynare structure describing the model
%   options_:     Dynare structure describing the options
%
% OUTPUTS
%   sim_array_:   [endo_nbr by periods by replic] array of simulations
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013-2023 Johannes Pfeifer
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
% For a copy of the GNU General Public License
% see <http://www.gnu.org/licenses/>.
if ver_less_than(dynare_version,'5.0')
    fid=fopen([M_.fname,'_simul'],'r');
else
    fid=fopen([M_.dname filesep 'Output' filesep M_.fname,'_simul'],'r');
end
if fid<3
    error(['Cannot open ', M_.fname,'_simul']);
end
simulations=fread(fid,[M_.endo_nbr*options_.periods options_.simul_replic],'float64');
fclose(fid);
sim_array=reshape(simulations,[M_.endo_nbr options_.periods options_.simul_replic]);
