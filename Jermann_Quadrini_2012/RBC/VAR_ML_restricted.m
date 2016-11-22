function [minus_loglik,Sigma_U]=VAR_ML_restricted(par_opt_vector,par_vector,par_index,Y,Z)
%function [minus_loglik,Sigma_U]=VAR_ML_restricted(par_opt_vector,par_vector,par_index,Y,Z)
% computes the log-likelihood of a restricted VAR
% Inputs:
%   par_opt_vector: [n_opt by 1] vector of parameter over which to optimize
%   par_vector:     [(trend_dummy+constant_dummy)*nvars+nvars^2*nlags+nvars*(nvars-1)/2 by 1] vector of VAR parameters
%   par_index:      [n_opt by 1] vector of indices of par_opt_vector in par_vector
%   Y:              [nvars by T] Matrix of observed variables 
%   Z:              [nvars*nlags by T] Matrix of regressors
%
% Outputs:
%   minus_loglik    [scalar]        negative of the log-likelihood function
%   Sigma_U         [nvars*nvars]   covariance matrix of the shocks
%
% Copyright (C) 2014-2016 Johannes Pfeifer
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

par_vector(par_index,:)=par_opt_vector;

if size(Y,1)>size(Y,2)
    error('Wrong dimensions of Y')
end
[nvars, T]=size(Y);
Y=Y(:);

%build AR matrices
betta=par_vector(1:end-(nvars*(nvars+1))/2,1); 
L=vec2cholmat(par_vector(end-(nvars*(nvars+1))/2+1:end));

Sigma_U=L*L';
loglik=-nvars*T/2*log(2*pi)-T/2*log(det(Sigma_U))-0.5*(Y-kron(Z',eye(nvars))*betta)'*kron(eye(T),eye(nvars)/Sigma_U)*(Y-kron(Z',eye(nvars))*betta);
minus_loglik=-loglik; %use of a minimizer