function outvalue=solve_calibration(par_vec,betta,delta,hiring_cost_share,x,alfa,N,xi,U,varphi)
% function outvalue=solve_calibration(par_vec,betta,delta,hiring_cost_share,x,alfa,N,xi,U,varphi)
% Solves equations (22),(48)-(52) of Gali (2010) for the labor disutility chi, unemployment
% weight parameter psi, G and C, given the hiring_cost_share-target
%
% Copyright (C) 2016 Johannes Pfeifer
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
psi=par_vec(1);
chi=par_vec(2);
G=par_vec(3);
C=par_vec(4);

outvalue(1)=N^(1-alfa)-(C+delta*N*G); % (48)
outvalue(2)=(1-betta*(1-delta))*G-xi*((1-alfa)*N^(-alfa)-chi*C*(N+psi*U)^varphi); % (49) with (52) plugged in
outvalue(3)=(1-x)*xi*psi*chi*C*(N+psi*U)^varphi-(1-xi)*G*x; % (50) with (52) plugged in
outvalue(4)=G/hiring_cost_share-(xi*chi*C*(N+psi*U)^varphi+(1-xi)*(1-alfa)*N^(-alfa)); %(22)
end