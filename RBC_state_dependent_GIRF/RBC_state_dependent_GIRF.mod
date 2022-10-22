/*
 * This file uses a basic RBC model with government spending and TFP shocks to
 * demonstrate how to compute Generalized Impulse Response Functions using 
 * Dynare's simult_-function. The model is solved up to second order to allow 
 * for non-linearities.
 * 
 * Two exercises are considered:
 *      1. A comparison of the absolute output response to a 1% of GDP shock that is either 
 *          positive or negative. The positive/negative shock is imposed at the ergodic mean.
 *          It shows that in nonlinear models, the sign of the shock matters.
 *      2. A comparison of the absolute output response to a positive 1% of GDP 
 *          shock at the ergodic mean and when capital is 10% below the steady state.
 *          It shows that in nonlinear models, the point in the state space where the shock happens 
 *          matters.
 *
* 
 * Notes:
 * - This implementation has been tested with Dynare 5.2
 *
 * This implementation was written by Johannes Pfeifer. 
 */

/*
 * Copyright (C) 2015-2022 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */

var y c k l z ghat r w invest; 
varexo eps_z eps_g;

parameters beta psi sigma delta alpha rho gammax rhog gshare l_ss k_ss i_ss y_ss g_ss c_ss n x k_y i_y;

sigma=5; %use high risk aversion to make model more nonlinear
alpha = 0.33;
i_y=0.25; %investment to output ration
k_y=10.4; %capital output ratio
x=0.0055; %average output growth
n=0.0027; %population growth
rho     = 0.97; %autocorrelation TFP
rhog    = 0.98; %autocorrelation government spending
gshare  = 0.2038; %share of government spending in GDP


model; 
  [name='FOC Labor (log utility)']
  psi*exp(c)^sigma*1/(1-exp(l)) = (1-alpha)*exp(z)*(exp(k(-1))/exp(l))^alpha;
  [name='Euler Equation']
  exp(c)^(-sigma) = beta/gammax*exp(c(+1))^(-sigma)*(alpha*exp(z(+1))*(exp(k)/exp(l(+1)))^(alpha-1)+(1-delta));
  [name='Resource constraint']
  gammax*exp(k) = exp(y)-exp(c)+(1-delta)*exp(k(-1))-g_ss*exp(ghat);
  [name='Production function']
  exp(y) = exp(z)*exp(k(-1))^alpha*exp(l)^(1-alpha);
  [name='LOM TFP']
  z = rho*z(-1)+eps_z;
  [name='LOM government spending']
  ghat = rhog*ghat(-1)+eps_g;
  [name='real wage']
  exp(w) = (1-alpha)*exp(y)/exp(l);
  [name='annualized interest rate']
  r =4* alpha*exp(y)/exp(k(-1)); // annualized
  [name='Definition investment']
  exp(invest)=exp(y)-exp(c)-g_ss*exp(ghat);
end;

%----------------------------------------------------------------
%  set steady state values and calibrate the model to steady state labor of 0.33, i.e. compute the corresponding steady state values
%  and the labor disutility parameter by hand
%---------------------------------------------------------------


steady_state_model;
    gammax=(1+n)*(1+x);
    delta=i_y/k_y-x-n-n*x;
    beta=(1+x)*(1+n)/(alpha/k_y+(1-delta));
    l_ss=0.33;
    k_ss = ((1/beta*(1+n)*(1+x)-(1-delta))/alpha)^(1/(alpha-1))*l_ss; 
    i_ss = (x+n+delta+n*x)*k_ss;
    y_ss=k_ss^alpha*l_ss^(1-alpha);
    g_ss=gshare*y_ss;
    c_ss = (1-gshare)*k_ss^(alpha)*l_ss^(1-alpha)-i_ss;
    psi=(1-alpha)*(k_ss/l_ss)^alpha*(1-l_ss)/c_ss^sigma;
    invest=log(i_ss);
    w = log((1-alpha)*y_ss/l_ss);
    r = 4*alpha*y_ss/k_ss;
    y = log(y_ss);
    k = log(k_ss);
    c = log(c_ss);
    l = log(l_ss);
    z = 0; 
    ghat =0;
end;

%----------------------------------------------------------------
%  set shock variances
%---------------------------------------------------------------

shocks;
    var eps_z=0.0068^2;
    var eps_g=0.0105^2;
end;

steady;
check;

%----------------------------------------------------------------
% compute policy function at second order
%----------------------------------------------------------------

stoch_simul(order = 2,irf=0);


%%%%%%%%%%%%%%%%%% Generate GIRF with positive 1% shock at ergodic mean %%%%

%define options
irf_periods=20; %IRF should have 20 periods
drop_periods=200; %drop 200 periods in simulation as burnin
irf_replication=1000; %take GIRF average over 1000 periods


impulse_vec=zeros(1,M_.exo_nbr); %initialize impulse vector to 0
impulse_vec(strmatch('eps_g',M_.exo_names,'exact'))=0.01/gshare; %impulse of 1 percent of GDP shock to g

starting_point=oo_.dr.ys; %define starting point of simulations; here: start at steady state and use 200 periods burnin to get to ergodic distribution

%% initialize shock matrices to 0
shocks_baseline = zeros(irf_periods+drop_periods,M_.exo_nbr); %baseline
shocks_impulse = shocks_baseline;

%%  eliminate shocks with 0 variance
i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 )); %finds shocks with 0 variance
nxs = length(i_exo_var); %number of those shocks
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));%get Cholesky of covariance matrix to generate random numbers

IRF_mat=NaN(M_.endo_nbr,irf_periods,irf_replication); %initialize matrix
for irf_iter = 1: irf_replication
    shocks_baseline(:,i_exo_var) = randn(irf_periods+drop_periods,nxs)*chol_S; %generate baseline shocks
    shocks_impulse = shocks_baseline; %use same shocks in impulse simulation
    shocks_impulse(drop_periods+1,:) = shocks_impulse(drop_periods+1,:)+impulse_vec; %add deterministic impulse
    y_baseline = simult_(M_,options_,starting_point,oo_.dr,shocks_baseline,options_.order); %baseline simulation
    y_shock = simult_(M_,options_,starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
    IRF_mat(:,:,irf_iter) = (y_shock(:,M_.maximum_lag+drop_periods+1:end)-y_baseline(:,M_.maximum_lag+drop_periods+1:end)); %add up difference between series
end
GIRF_positive=mean(IRF_mat,3); %take average


%%%%%%%%%%%%%%%%%% Generate GIRF with negative 1% shock at ergodic mean %%%%


impulse_vec=zeros(1,M_.exo_nbr); %initialize impulse vector to 0
impulse_vec(strmatch('eps_g',M_.exo_names,'exact'))=-0.01/gshare; %impulse of minus 1 percent of GDP shock to g

starting_point=oo_.dr.ys; %define starting point of simulations; here: start at steady state and use 200 periods burnin to get to ergodic distribution

%% initialize shock matrices to 0
shocks_baseline = zeros(irf_periods+drop_periods,M_.exo_nbr); %baseline
shocks_impulse = shocks_baseline;

%%  eliminate shocks with 0 variance
i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 )); %finds shocks with 0 variance
nxs = length(i_exo_var); %number of those shocks
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));%get Cholesky of covariance matrix to generate random numbers

IRF_mat=NaN(M_.endo_nbr,irf_periods,irf_replication);
for irf_iter = 1: irf_replication
    shocks_baseline(:,i_exo_var) = randn(irf_periods+drop_periods,nxs)*chol_S; %generate baseline shocks
    shocks_impulse = shocks_baseline; %use same shocks in impulse simulation
    shocks_impulse(drop_periods+1,:) = shocks_impulse(drop_periods+1,:)+impulse_vec; %add deterministic impulse
    y_baseline = simult_(M_,options_,starting_point,oo_.dr,shocks_baseline,options_.order); %baseline simulation
    y_shock = simult_(M_,options_,starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
    IRF_mat(:,:,irf_iter) = (y_shock(:,M_.maximum_lag+drop_periods+1:end)-y_baseline(:,M_.maximum_lag+drop_periods+1:end)); %add up difference between series
end
GIRF_negative=mean(IRF_mat,3); %take average


figure('Name','GIRFs to positive and negative G shock at ergodic mean')
subplot(2,1,1)
plot(1:irf_periods,GIRF_positive(strmatch('ghat',M_.endo_names,'exact'),:),'b-',1:irf_periods,-GIRF_negative(strmatch('ghat',M_.endo_names,'exact'),:),'r--')
title('G')
subplot(2,1,2)
plot(1:irf_periods,GIRF_positive(strmatch('y',M_.endo_names,'exact'),:),'b-',1:irf_periods,-GIRF_negative(strmatch('y',M_.endo_names,'exact'),:),'r--')
title('|Y|')
legend('Positive G-Shock','Minus Response to Negative G-shock')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate GIRF with positive 1% shock at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% with capital 10% below steady state

%define options
irf_periods=20; %IRF should have 20 periods
drop_periods=0; %drop 0 periods in simulation as burnin
irf_replication=1000; %take GIRF average over 1000 periods

impulse_vec=zeros(1,M_.exo_nbr); %initialize impulse vector to 0
impulse_vec(strmatch('eps_g',M_.exo_names,'exact'))=0.01/gshare; %impulse of 1 percent of GDP shock to g


starting_point=oo_.dr.ys; %define starting point of simulations; here: start at steady state
starting_point(strmatch('k',M_.endo_names,'exact'))=0.9*starting_point(strmatch('k',M_.endo_names,'exact')); %set capital to 10% below steady state

%% initialize shock matrices to 0
shocks_baseline = zeros(irf_periods+drop_periods,M_.exo_nbr); %baseline
shocks_impulse = shocks_baseline;

%%  eliminate shocks with 0 variance
i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 )); %finds shocks with 0 variance
nxs = length(i_exo_var); %number of those shocks
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));%get Cholesky of covariance matrix to generate random numbers

IRF_mat=NaN(M_.endo_nbr,irf_periods,irf_replication);
for irf_iter = 1: irf_replication
    shocks_baseline(:,i_exo_var) = randn(irf_periods+drop_periods,nxs)*chol_S; %generate baseline shocks
    shocks_impulse = shocks_baseline; %use same shocks in impulse simulation
    shocks_impulse(drop_periods+1,:) = shocks_impulse(drop_periods+1,:)+impulse_vec; %add deterministic impulse
    y_baseline = simult_(M_,options_,starting_point,oo_.dr,shocks_baseline,options_.order); %baseline simulation
    y_shock = simult_(M_,options_,starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
    IRF_mat(:,:,irf_iter) = (y_shock(:,M_.maximum_lag+drop_periods+1:end)-y_baseline(:,M_.maximum_lag+drop_periods+1:end)); %add up difference between series
end
GIRF_positive_low_capital=mean(IRF_mat,3); %take average


figure('Name','GIRFs to positive G shock at different points')
subplot(2,1,1)
plot(1:irf_periods,GIRF_positive(strmatch('ghat',M_.endo_names,'exact'),:),'b-',1:irf_periods,GIRF_positive_low_capital(strmatch('ghat',M_.endo_names,'exact'),:),'r--')
title('G')
subplot(2,1,2)
plot(1:irf_periods,GIRF_positive(strmatch('y',M_.endo_names,'exact'),:),'b-',1:irf_periods,GIRF_positive_low_capital(strmatch('y',M_.endo_names,'exact'),:),'r--')
title('|Y|')
legend('Positive G-Shock at ergodic mean','Positive G-Shock with capital 10% below SS')