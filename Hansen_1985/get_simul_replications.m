function [sim_array]=get_simul_replications(DynareModel,options_)
% function [sim_array]=get_simul_replications(DynareModel,options_)
% reads the simulation replications into a three-dimensional array with
% endogenous variables along the first, simulation periods along the
% second, and replications along the third dimension
%
% INPUTS
%   DynareModel:  Dynare structure describing the model
%   options_:     Dynare structure describing the options
%
% OUTPUTS
%   sim_array_:   [endo_nbr by periods by replic] array of simulations
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013-15 Johannes Pfeifer
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

fid=fopen([DynareModel.fname,'_simul'],'r');
if fid<3
    error(['Cannot open ', DynareModel.fname,'_simul']);
end
simulations=fread(fid,[DynareModel.endo_nbr*options_.periods options_.simul_replic],'float64');
fclose(fid);
sim_array=reshape(simulations,[DynareModel.endo_nbr options_.periods options_.simul_replic]);
