function [F, Q, H]=build_companionform(A,p,Q_start,H_start)
%[F, Q, H]=build_companionform(A_start,numberoflags,Q_start,H_start)
% Builds companion form to generate AR-1 process from AR-p for state-space
% system
%
% Inputs:
%   - A             [nstates by nstates*numberoflags] matrix containing the p AR-matrices
%   - p             [scalar]    number of lags
%   - Q_start       [nstates by nstates] covariance matrix of structural shocks
%   - H_start       [nobservables by nstates] Matrix mapping states
%                   into observables

% Outputs:
%   - F             [numberoflags*nstates by nstates*numberoflags] companion form matrix of
%                       the state transition equation
%   - Q             [numberoflags*nstates by numberoflags*nstates] covariance matrix of
%                       structural shocks in companion form
%   - H             [nobservables by nstates*nobservables] Matrix mapping states
%                   into observables

% Copyright (C) 2015-2017 Johannes Pfeifer
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


nstates=size(A,1);
if size(A,2)~=nstates*p
    error('Matrix size not conformable with lags');
end
F=zeros(nstates*p,nstates*p);
for lag_iter=1:p 
    F(1:nstates,(lag_iter-1)*nstates+1:lag_iter*nstates)=A(:,(lag_iter-1)*nstates+1:lag_iter*nstates);
end
F(nstates+1:nstates*p,1:nstates*(p-1))=eye(nstates*(p-1),nstates*(p-1));

if nargin>=3
    Q=zeros(nstates*p,nstates*p);
    Q(1:nstates,1:nstates)=Q_start;
end
if nargin>=4
    nobservables=size(H_start,1);
    H=zeros(size(H_start,1),nstates*p);
    H(1:nobservables,1:nstates)=H_start;
end
