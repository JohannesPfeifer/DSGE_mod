function y = simult_nonsymmetric_order3(M_,oo_,options_,exo)
% y = simult_nonsymmetric_order3(M_,oo_,options_)
% -------------------------------------------------------------------------
% This function performs a stochastic simulation for variables declared as
% VAROBS using a third-order (unpruned) perturbation approximation taking
% into account that the shocks stem from a non-Gaussian distribution
% specified by nonzero third-order product moments as this makes ghs3 nonzero.
% If the second-order product moments also change, one needs to re-compute
% the perturbation matrices ghs2, ghs3, ghxss, and ghuss as well.
% See for example perturbation_solver_nonsymmetric_order3.m for such a function.
% -------------------------------------------------------------------------
% NOTATION
% - x denote PREVIOUS state variables (predetermined and mixed)
% - y denote CURRENT VAROBS variables
% - u denote CURRENT shocks
% - s denotes perturbation parameter
% =========================================================================
% INPUTS
%  * M_:                     [structure] information about model
%  * oo_:                    [structure] storage for results (particularly lower-order perturbation matrices)
%  * exo:                    [options_.periods by M_.exo_nbr] matrix of random shocks
% -------------------------------------------------------------------------
% OUTPUTS
%  * y:                      [length(options_.varobs) by options_.periods] simulated data for VAROBS variables
% =========================================================================
% Copyright (C) 2022 Willi Mutschler (@wmutschl, willi@mutschler.eu)
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
% see <https://www.gnu.org/licenses/>.
% =========================================================================

% create joint index for state variables x and varobs variables y (we will focus only on VAROBS variables below)
indx = M_.nstatic + (1:M_.nspred); % state variables in DR order
indy = [];
for j=1:length(options_.varobs) %use for-loop to keep same ordering as declared in varobs
    indy = [indy find(ismember(M_.endo_names(oo_.dr.order_var),options_.varobs{j}))]; % reporting variables in DR order
end
indxy = [indx(:);indy(:)]; % joint index for state and observable variables
xy_ss = oo_.dr.ys(oo_.dr.order_var); xy_ss = xy_ss(indxy); % steady-state of x and y in DR order

% select rows in perturbation matrices for state and observable variables only
gx = oo_.dr.ghx(indxy,:);
gu = oo_.dr.ghu(indxy,:);
gxx = oo_.dr.ghxx(indxy,:);
gxu = oo_.dr.ghxu(indxy,:);
guu = oo_.dr.ghuu(indxy,:);
gss = oo_.dr.ghs2(indxy,:);
gxxx = oo_.dr.ghxxx(indxy,:);
gxxu = oo_.dr.ghxxu(indxy,:);
gxuu = oo_.dr.ghxuu(indxy,:);
gxss = oo_.dr.ghxss(indxy,:);
guuu = oo_.dr.ghuuu(indxy,:);
guss = oo_.dr.ghuss(indxy,:);
gsss = oo_.dr.ghs3(indxy,:);
    
% initialize
xhat = zeros(length(indx),1);
y = zeros(length(indy),options_.periods);
% do the simulations (one can also replace kron() with Dynare's mex A_times_B_kronecker_C which computes A*kron(B,C) more efficiently)
%hh = dyn_waitbar(0,'Simulating the model');
%set(hh,'Name','Simulating the model');
for t = 1:options_.periods
    %if mod(t,100000)==1
    %    dyn_waitbar(t/options_.periods, hh, sprintf('Period %u of %u', t,options_.periods));
    %end
    u = exo(t,:)';
    xy = xy_ss ...
       + gx*xhat + gu*u ...
       + 1/2*gxx*kron(xhat,xhat) + gxu*kron(xhat,u) + 1/2*guu*kron(u,u) + 1/2*gss ...
       + 1/6*gxxx*kron(xhat,kron(xhat,xhat)) ...
       + 1/6*guuu*kron(u,kron(u,u)) ...
       + 3/6*gxxu*kron(xhat,kron(xhat,u)) ...
       + 3/6*gxuu*kron(xhat,kron(u,u)) ...
       + 3/6*gxss*xhat ...
       + 3/6*guss*u ...
       + 1/6*gsss ...
       ;
    xhat = xy(1:length(indx)) - xy_ss(1:length(indx)); % update states (in deviation from steady-state)
    y(:,t) = xy((length(indx)+1):end); % update observables
end
%dyn_waitbar_close(hh);