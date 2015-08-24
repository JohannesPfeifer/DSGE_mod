function [x,y,policyfun]=plot_policy_fun(state_name,state_range,plot_var_name,y0,plot_dummy)
% Plots the policy function
%
% INPUTS
%    state_name     [string]   string denoting the state in whose dependence to plot the
%                              policy function (x-axis)
%    state_range    [double]   k*1 vector, providing the numerical values for
%                              the state grid
%    plot_var_name  [string]   string denoting the name of the variable for which to plot the
%                              policy function
%    y0             [double]   n*1 vector, initial value, defaults to SS, if not specified (n is the number of declared endogenous variables plus the number 
%                              of auxilliary variables for lags and leads)
%    plot_dummy     [scalar]   dummy whether plot should be created
%
% OUTPUTS
%    x               [double]   1*k vector of numerical values for the state grid.
%    y               [double]   1*k vector of numerical values for the endogenous variable to be plotted.
%    policyfun       [double]   n*k matrix of numerical values for all endogenous variable
% SPECIAL REQUIREMENTS
%    none

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

global M_ options_ oo_

if nargin<4 || isempty(y0)
    y0=oo_.dr.ys;
end
if nargin<5
    plot_dummy=1; 
end
n_points=length(state_range);

policyfun=zeros(M_.endo_nbr,n_points);

if options_.block
    if M_.maximum_lag > 0
        k2 = oo_.dr.state_var;
    else
        k2 = [];
    end;
    order_var = 1:M_.endo_nbr;
    oo_.dr.order_var = order_var;
else
    k2 = oo_.dr.kstate(find(oo_.dr.kstate(:,2) <= M_.maximum_lag+1),[1 2]);
    k2 = k2(:,1)+(M_.maximum_lag+1-k2(:,2))*M_.endo_nbr;
    order_var = oo_.dr.order_var;
end;

if ismember(state_name,M_.endo_names(order_var(k2),:))
    xlabelstring=[state_name,'(-1)'];
elseif ismember(state_name,M_.exo_names)
    xlabelstring=[state_name];
elseif ismember(state_name,M_.endo_names) && ~ismember(state_name,M_.endo_names(order_var(k2),:))
    error([state_name,' is not a state variable'])
else
    error(['Unknown Variable ',state_name])
end
if ~ismember(plot_var_name,M_.endo_names)
    error(['Unknown Variable ',plot_var_name,' to plot'])    
end

for ii=1:n_points
    y0_temp=y0;
    shock_mat=zeros(1,M_.exo_nbr);
    if ismember(state_name,M_.endo_names(order_var(k2),:))
        y0_temp(strmatch(state_name,M_.endo_names,'exact'),1)=state_range(1,ii);
    elseif ismember(state_name,M_.exo_names)
        shock_mat(1,strmatch(state_name,M_.exo_names,'exact'),1)=state_range(1,ii);
    end
    temp=simult_(y0_temp,oo_.dr,shock_mat,options_.order);
    policyfun(:,ii)=temp(:,2);
end
x=state_range;
y=policyfun(strmatch(plot_var_name,M_.endo_names,'exact'),:);
if plot_dummy
    figure('Name',['Policy function for ',plot_var_name])
    plot(state_range,y);
    xlabel(xlabelstring,'Interpreter','none')
    ylabel(plot_var_name,'Interpreter','none')
end

