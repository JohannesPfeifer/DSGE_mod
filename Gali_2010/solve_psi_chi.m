function outvalue=solve_psi_chi(par_vec,betta,delta,Gamma,x,gammma,alfa,N,C,xi,U,varphi)
% function outvalue=solve_psi_chi(par_vec,betta,delta,Gamma,x,gammma,alfa,N,C,xi,U,varphi)
% Solves equations (49)-(50) of Gali (2010) for the labor disutility chi and unemployment
% weight parameter psi
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
outvalue(1)=(1-betta*(1-delta))*Gamma*x^gammma-xi*((1-alfa)*N^(-alfa)-chi*C*(N+psi*U)^varphi); % (49) with (52) plugged in
outvalue(2)=(1-x)*xi*psi*chi*C*(N+psi*U)^varphi-(1-xi)*Gamma*x^(1+gammma); % (50) with (52) plugged in
end