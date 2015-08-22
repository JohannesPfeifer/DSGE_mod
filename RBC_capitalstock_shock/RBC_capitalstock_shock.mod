/*
 * This file implements a simple RBC model with a time t shock to the capital stock.
 * Basically, the capital stock used in time t is not completey predetemined anymore. 
 * While the model conforms to Dynare's end of period stock notation, the production 
 * function uses k instead of k(-1), because the effective capital stock used at time t 
 * for production can change at time to due to a shock.
 * Similarly, the timing for k in all other equations changes to being contemporaneous. 
 * At the same time, the law of motion for capital has to be adjusted accordingly to  
 * reflect the presence of the shock. The capital stock at time t is still
 * determined by yesterday's capital and yesterday's investment, i.e. it is 
 * predetermined in terms of endogenous variables, but it is also affected by a time t 
 * exogenous shock named eps_cap that leads to a percentage drop in capital:
 *
 * k = exp(-eps_cap)*(invest(-1)+(1-delta)*k(-1));
 * 
 * The model uses the variable redefinition x=log(X) <=> exp(x)=X to use Dynare's
 * linearization capabilities to perform a log-linearization. Thus, all impulse 
 * response functions are in percentage deviations from the deterministic steady state.
 * 
 * This implementation was written by Johannes Pfeifer. In
 * case you spot mistakes, email us at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.

 */

/*
 * Copyright (C) 2013-15 Johannes Pfeifer
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

%----------------------------------------------------------------
% define variables 
%----------------------------------------------------------------

var y c k l z invest;
varexo eps_z eps_cap;

%----------------------------------------------------------------
% define parameters
%----------------------------------------------------------------

parameters beta psi delta alpha rho i_y k_y l_ss k_ss i_ss y_ss c_ss;

%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------
alpha   = 0.33;
i_y = 0.25;
k_y = 10.4;
rho     = 0.97;  

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model; 
  psi*exp(c)/(1-exp(l)) = (1-alpha)*exp(z)*(exp(k)/exp(l))^alpha;
  1/exp(c)= beta/exp(c(+1))*(alpha*exp(z(+1))*(exp(k(+1))/exp(l(+1)))^(alpha-1)+(1-delta));
  exp(k) = exp(-eps_cap)*(exp(invest(-1))+(1-delta)*exp(k(-1)));
  exp(y) = exp(z)*exp(k)^alpha*exp(l)^(1-alpha);
  z = rho*z(-1)+eps_z;
  exp(invest)=exp(y)-exp(c);
end;

%----------------------------------------------------------------
%  set steady state values
%---------------------------------------------------------------


steady_state_model;
    //Do Calibration
    delta = i_y/k_y;
    beta = 1/(alpha/k_y+(1-delta));
    l_ss=0.33;
    k_ss = ((1/beta-(1-delta))/alpha)^(1/(alpha-1))*l_ss; 
    i_ss = delta*k_ss;
    y_ss=k_ss^alpha*l_ss^(1-alpha);
    c_ss = k_ss^(alpha)*l_ss^(1-alpha)-i_ss;
    psi=(1-alpha)*(k_ss/l_ss)^alpha*(1-l_ss)/c_ss;
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
var eps_z = 1;
var eps_cap = 1; 
end;

resid(1);

steady;

check;

%----------------------------------------------------------------
% generate IRFs
%----------------------------------------------------------------

stoch_simul(order = 1,irf=20);


