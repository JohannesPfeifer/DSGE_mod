/*
 * This file presents a baseline RBC model with TFP shocks where a distortionary labor tax rate 
 * is levied and then rebated to households via lump-sum taxation
 * The mod-file then shows how to maximize unconditional welfare (mean welfare) by choosing the optimal
 * labor tax rate, which - unsurprisingly - is 0. It does so by running an optimizer over a function 
 * return_welfare.m setting the parameter to be optimized and returning the objective value.
 *
 * The model setup is similar to the one described in Handout_RBC_model.pdf and 
 *  resembles the one in King/Rebelo (1999): 
 *  Resuscitating Real Business Cycles, Handbook of Macroeconomics, Volume 1 and
 *  Romer (2012), Advanced macroeconomics, 4th edition.
 * 
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2017 Johannes Pfeifer
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

//****************************************************************************
//Define variables
//****************************************************************************

var y           ${y}$ (long_name='output')
    c           ${c}$ (long_name='consumption')
    k           ${k}$ (long_name='capital')
    l           ${l}$ (long_name='hours')
    z           ${z}$ (long_name='TFP')
    r           ${r}$ (long_name='annualized interest rate')
    w           ${w}$ (long_name='real wage')
    invest      ${i}$ (long_name='investment') 
    log_y       ${\log(y)}$ (long_name='log output')
    log_k       ${\log(k)}$ (long_name='log capital stock')
    log_c       ${\log(c)}$ (long_name='log consumption')
    log_l       ${\log(l)}$ (long_name='log labor')
    log_w       ${\log(w)}$ (long_name='log real wage')
    log_invest  ${\log(i)}$ (long_name='log investment')
    W           ${W}$ (long_name='Welfare')
    ;

varexo eps_z ${\varepsilon_z}$ (long_name='TFP shock');
    
parameters 
    betta   ${\beta}$   (long_name='discount factor')
    psii    ${\psi}$    (long_name='labor disutility parameter')
    siggma  ${\sigma}$  (long_name='risk aversion')
    delta   ${\delta}$  (long_name='depreciation rate')
    alppha  ${\alpha}$  (long_name='capital share')
    rhoz    ${\rho_z}$  (long_name='persistence TFP shock')
    gammax  ${\gamma_x}$ (long_name='composite growth rate')
    tau_n   ${\tau^n}$  (long_name='labor tax rate')
    n       ${n}$       (long_name='population growth')
    x       ${x}$       (long_name='technology growth (per capita output growth)')
    i_y     ${\frac{I}{Y}}$ (long_name='investment-output ratio')
    k_y     ${\frac{K}{Y}}$ (long_name='capital-output ratio')
    ;

//****************************************************************************
//Set parameter values
//****************************************************************************

siggma=2;                % risk aversion
alppha= 0.33;            % capital share
i_y=0.25;               % investment-output ration
k_y=10.4;               % capital-output ratio
x=0.0055;               % technology growth (per capita output growth)
n=0.0027;               % population growth
rhoz=0.97;              %technology autocorrelation base on linearly detrended Solow residual
tau_n=0.2;
psii=3.4880;

//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************

model;
[name='Euler equation']
c^(-siggma)=betta/gammax*c(+1)^(-siggma)*
    (alppha*exp(z(+1))*(k/l(+1))^(alppha-1)+(1-delta));
[name='Recursive Welfare definition']    
W=c^(1-siggma)/(1-siggma)+psii*log(1-l)+betta*W(+1);    
[name='Labor FOC']
psii*c^siggma*1/(1-l)=(1-tau_n)*w;
[name='Law of motion capital'] 
gammax*k=(1-delta)*k(-1)+invest;
[name='resource constraint']
y=invest+c;
[name='production function']
y=exp(z)*k(-1)^alppha*l^(1-alppha);
[name='real wage/firm FOC labor']
w=(1-alppha)*y/l;
[name='annualized real interest rate/firm FOC capital']
r=4*alppha*y/k(-1);
[name='exogenous TFP process']
z=rhoz*z(-1)+eps_z;
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

//****************************************************************************
//set shock variances
//****************************************************************************

shocks;
    var eps_z=0.66^2;
end;

//****************************************************************************
//check the starting values for the steady state
//****************************************************************************

resid;

//****************************************************************************
// compute steady state given the starting values
//****************************************************************************

steady;
//****************************************************************************
// check Blanchard-Kahn-conditions
//****************************************************************************

check;

//****************************************************************************
// compute policy function at second order
//****************************************************************************

stoch_simul(order=2,irf=0) W;

options_.noprint=1; %shut off error messages and output
[tau_hat,W_hat] = fmincon(@return_welfare,0.2,[],[],[],[],0,1);
options_.noprint=0; 

%compute output under optimal policy
set_param_value('tau_n',tau_hat); %set tax rate to optimal one
stoch_simul(order=2,irf=0) W c l;