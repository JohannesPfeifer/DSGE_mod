/*
 * This file presents a baseline RBC model with government spending shocks where the persistence of the 
 *  government spending shock is estimated via impulse response function (IRF) matching. 
 *
 * Notes:
 *  - The empirical IRFs are estimated using the Blanchard/Perotti (2002) approach. Of course the RBC
 *      model is not capable of generating the consumption increase after a government spending shock. For that
 *      reason, this mod-file only targets the IRFs for G and Y.
 *  - The weighting matrix uses a diagonal matrix with the inverse of the pointwise IRF variances on the main 
 *      diagonal. The same approach has for example been used in Christiano/Eichenbaum/Evans (2005). To estimate 
 *      the variances, a simple residual bootstrap is performed
 *  - The empirical IRFs and model IRFs use an impulse size of 1 percent. Thus, there is no uncertainty about the 
 *      initial impact. The IRF matching therefore only targets the G-response starting in the second period.
 *  - Note that for the current model, the number of IRFs exceeds the number of VAR parameters. Therefore,
 *      the distribution of the estimator will be non-standard, see Guerron-Quintana/Inoue/Kilian (2016), 
 *      http://dx.doi.org/10.1016/j.jeconom.2016.09.009
 *  - The mod-file also shows how to estimate an AR(2)-process by working with the roots of the autoregressive
 *      process instead of the coefficients. This allows for easily restricting the process to the stability region and 
 *      would allow specifying e.g. a beta prior for both roots as was done in Born/Peter/Pfeifer (2013), Fiscal news 
 *      and macroeconomic volatility, https://doi.org/10.1016/j.jedc.2013.06.011
 *  - The dataset was donwloaded using the FRED plugin and can be easily updated with it.    
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2016-17 Johannes Pfeifer
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
@#define IRF_periods=80
@#define CMAES=0

var y           ${y}$ (long_name='output')
    c           ${c}$ (long_name='consumption')
    k           ${k}$ (long_name='capital')
    l           ${l}$ (long_name='hours')
    z           ${z}$ (long_name='TFP')
    ghat        ${\hat g}$ (long_name='government spending')
    r           ${r}$ (long_name='annualized interest rate')
    w           ${w}$ (long_name='real wage')
    invest      ${i}$ (long_name='investment') 
    log_y       ${\log(y)}$ (long_name='log output')
    log_k       ${\log(k)}$ (long_name='log capital stock')
    log_c       ${\log(c)}$ (long_name='log consumption')
    log_l       ${\log(l)}$ (long_name='log labor')
    log_w       ${\log(w)}$ (long_name='log real wage')
    log_invest  ${\log(i)}$ (long_name='log investment')
    ;

varexo eps_z ${\varepsilon_z}$ (long_name='TFP shock')
       eps_g ${\varepsilon_g}$ (long_name='government spending shock')
    ;

%----------------------------------------------------------------
% define parameters
%----------------------------------------------------------------

parameters 
    beta    ${\beta}$   (long_name='discount factor')
    psi     ${\psi}$    (long_name='labor disutility parameter')
    sigma   ${\sigma}$  (long_name='risk aversion')
    delta   ${\delta}$  (long_name='depreciation rate')
    alpha   ${\alpha}$  (long_name='capital share')
    rhoz    ${\rho_z}$  (long_name='persistence TFP shock')
    root_g_1    ${\rho_g}$  (long_name='persistence G shock')
    root_g_2    ${\rho_g}$  (long_name='persistence G shock')
    gammax  ${\gamma_x}$ (long_name='composite growth rate')
    gshare  ${\frac{G}{Y}}$ (long_name='government spending share')
    n       ${n}$       (long_name='population growth')
    x       ${x}$       (long_name='technology growth (per capita output growth)')
    i_y     ${\frac{I}{Y}}$ (long_name='investment-output ratio')
    k_y     ${\frac{K}{Y}}$ (long_name='capital-output ratio')
    g_ss    ${\bar G}$ (long_name='government spending in steady state')
    ;

%----------------------------------------------------------------
% set parameter values 
%----------------------------------------------------------------
sigma=1;
alpha= 0.33;
i_y=0.25;
k_y=10.4;
x=0.0055;
n=0.0027;
rhoz=0.97;
root_g_1=0.9602;
root_g_2=0;
gshare=0.2038;

%----------------------------------------------------------------
% enter model equations
%----------------------------------------------------------------

model;
# rho_g_1= (root_g_1+root_g_2);
# rho_g_2= - root_g_1*root_g_2;
[name='Euler equation']
c^(-sigma)=beta/gammax*c(+1)^(-sigma)*
    (alpha*exp(z(+1))*(k/l(+1))^(alpha-1)+(1-delta));
[name='Labor FOC']
psi*c^sigma*1/(1-l)=w;
[name='Law of motion capital'] 
gammax*k=(1-delta)*k(-1)+invest;
[name='resource constraint']
y=invest+c+g_ss*exp(ghat);
[name='production function']
y=exp(z)*k(-1)^alpha*l^(1-alpha);
[name='real wage/firm FOC labor']
w=(1-alpha)*y/l;
[name='annualized real interst rate/firm FOC capital']
r=4*alpha*y/k(-1);
[name='exogenous TFP process']
z=rhoz*z(-1)+eps_z;
[name='government spending process']
ghat=rho_g_1*ghat(-1)+rho_g_2*ghat(-2)+eps_g;
[name='Definition log output']
log_y = log(y);
[name='Definition log capital']
log_k = log(k);
[name='Definition log consumption']
log_c = log(c);
[name='Definition log hours']
log_l = log(l);
[name='Definition log wage']
log_w = log(w);
[name='Definition log investment']
log_invest = log(invest);
end;

%----------------------------------------------------------------
%  set steady state values
%---------------------------------------------------------------

steady_state_model;
    //Do Calibration
    //calibrate the model to steady state labor of 0.33,i.e. compute the corresponding steady state values
    // and the labor disutility parameter by hand;
    gammax=(1+n)*(1+x);
    delta=i_y/k_y-x-n-n*x;
    beta=(1+x)*(1+n)/(alpha/k_y+(1-delta));
    l=0.33;
    k = ((1/beta*(1+n)*(1+x)-(1-delta))/alpha)^(1/(alpha-1))*l; 
    invest = (x+n+delta+n*x)*k;
    y=k^alpha*l^(1-alpha);
    g=gshare*y;
    g_ss=g;
    c = (1-gshare)*k^(alpha)*l^(1-alpha)-invest;
    psi=(1-alpha)*(k/l)^alpha*(1-l)/c^sigma;
    w = (1-alpha)*y/l;
    r = 4*alpha*y/k;
    log_y = log(y);
    log_k = log(k);
    log_c = log(c);
    log_l = log(l);
    log_w = log(w);
    log_invest = log(invest);
    z = 0; 
    ghat =0;
end;

%----------------------------------------------------------------
%  set shock variances
%---------------------------------------------------------------

shocks;
// var eps_z = 0.0066^2;
var eps_g = 1; 
end;

resid(1);

steady;

check;

%----------------------------------------------------------------
% generate IRFs and compute model moments
%----------------------------------------------------------------
stoch_simul(order = 1,irf=@{IRF_periods}) log_y log_c ghat;


%% get empirical IRFs and weighting matrix
[IRF_empirical,IRF_weighting,IRF_quantiles]=get_empirical_IRFs(@{IRF_periods})
         

x_start=[root_g_1 root_g_2]; %use calibration as starting point
%make sure Dynare does not print out stuff during runs
options_.nomoments=1;
options_.nofunctions=1;
options_.nograph=1;
options_.verbosity=0;

%set noprint option to suppress error messages within optimizer
options_.noprint=1;

@#if CMAES==0
    % set csminwel options
    H0 = 1e-2*eye(length(x_start)); %Initial Hessian 
    crit = 1e-8; %Tolerance
    nit = 1000;  %Number of iterations

    [fhat,x_opt_hat] = csminwel(@IRF_matching_objective,x_start,H0,[],crit,nit,IRF_empirical,IRF_weighting);
@#else
    %set CMAES options
    H0=0.2*ones(size(x_start,1),1)
    cmaesOptions = options_.cmaes;
    cmaesOptions.LBounds = [-1;-1];
    cmaesOptions.UBounds = [1;1];
    [x_opt_hat, fhat, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes('IRF_matching_objective',x_start,H0,cmaesOptions,IRF_empirical,IRF_weighting);
    x_opt_hat=BESTEVER.x;
@#endif

%get IRFs at the optimum and plot them
[fval, IRF_model]=IRF_matching_objective(x_opt_hat,IRF_empirical,IRF_weighting);

figure
subplot(2,1,1)
plot(1:options_.irf,IRF_empirical(:,1),1:options_.irf,IRF_model(:,1),1:options_.irf,IRF_quantiles(1,:,1),'r--',1:options_.irf,IRF_quantiles(1,:,2),'r--');
title('G')
subplot(2,1,2)
plot(1:options_.irf,IRF_empirical(:,2),1:options_.irf,IRF_model(:,2),1:options_.irf,IRF_quantiles(2,:,1),'r--',1:options_.irf,IRF_quantiles(2,:,2),'r--');
title('Y')
legend('Empirical','Model')


