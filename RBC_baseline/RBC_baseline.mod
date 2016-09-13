/*
 * This file presents a baseline RBC model with TFP and government spending shocks, calibrated to US data from
 *  1947Q4:2016Q1. The model setup is described in Handout_RBC_model.pdf and resembles the one in King/Rebelo (1999): 
 *  Resuscitating Real Business Cycles, Handbook of Macroeconomics, Volume 1 and
 *  Romer (2012), Advanced macroeconomics, 4th edition
 *  The driving processes are estimated as AR(1)-processes on linearly detrended data.
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2016 Johannes Pfeifer
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

var y c k l z ghat r w  invest log_y log_k log_c log_l log_w log_invest;
varexo eps_z eps_g;

//****************************************************************************
//Define parameters
//****************************************************************************

parameters beta psi sigma delta alpha rhoz rhog gammax
    gshare n x i_y k_y g_ss;

//****************************************************************************
//Set parameter values
//****************************************************************************

sigma=1;                // risk aversion
alpha= 0.33;            // capital share
i_y=0.25;               // investment-output ration
k_y=10.4;               // capital-output ratio
x=0.0055;               // technology growth (per capita output growth)
n=0.0027;               // population growth
rhoz=0.97;              //technology autocorrelation base on linearly detrended Solow residual
rhog=0.989;
gshare=0.2038;          //government spending share

//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************

model;
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
ghat=rhog*ghat(-1)+eps_g;
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
// Provide steady state values and calibrate the model to steady state labor of 0.33, 
// i.e. compute the corresponding steady state values
// and the labor disutility parameter by hand;
//****************************************************************************

steady_state_model;
    //Do Calibration
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

//****************************************************************************
//set shock variances
//****************************************************************************

shocks;
    var eps_z=0.66^2;
    var eps_g=1.04^2;
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
// compute policy function at first order, do IRFs and compute moments with HP-filter
//****************************************************************************

stoch_simul(order=1,irf=40,hp_filter=1600) log_y log_k log_c log_l log_w r z ghat;