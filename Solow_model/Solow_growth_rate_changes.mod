/*
 * This file uses Dynare's perfect foresight solver to study the transition 
 * behavior of simple a Solow-Swan economy with Cobb-Douglas production function
 * after unanticipated changes in technology or population growth. It starts the simulation 
 * in steady state for one growth rate and then considers two unanticipated
 * change at t=0:
 *  i)  a drop in the rate of technology growth from a positive value to 0
 *  ii) a drop in the rate of population growth from a positive value to 0
 * 
 * Notes:
 *  - Because the model is purely backward-looking, only an initval-block is
 *      required for the endogenous dynamics. 
 *  - The endval block's purpose here is 
 *      i)  to permanently change the growth rates (g or n) to a new value after
 *          the initial period 
 *      ii) to provide starting values for the perfect foresight solver
 *
 * This implementation was written by Johannes Pfeifer. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2014-2016 Johannes Pfeifer
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
var c           ${c}$ (long_name='consumption (intensive form)')
    k           ${k}$ (long_name='capital (intensive form)')
    y           ${y}$ (long_name='output (intensive form)')
    invest      ${i}$ (long_name='investment (intensive form)')
    log_c       ${\log c}$ (long_name='log consumption (intensive form)')
    log_k       ${\log k}$ (long_name='log capital (intensive form)')
    log_y       ${\log y}$ (long_name='log output (intensive form)')
    log_invest  ${\log i}$ (long_name='log investment (intensive form)')
    g_k_aggregate       ${\Delta K}$ (long_name='Growth rate aggregate capital')
    g_k_per_capita      ${\Delta \tilde k}$ (long_name='Growth rate capital per capita')
    g_k_intensive       ${\Delta k}$ (long_name='Growth rate capital in intensive form')
    ;

//****************************************************************************
//Define predetermined variables
//****************************************************************************
predetermined_variables k;

//****************************************************************************
//Define exogenous variables (n and g are time-varying in our experiment)
//****************************************************************************
varexo n    ${n}$ (long_name='population growth rate')
       g    ${g}$ (long_name='technology growth rate')
    ;

//****************************************************************************
//Define parameters
//****************************************************************************
parameters s    ${s}$ (long_name='saving rate')
    alpha       ${\alpha}$ (long_name='capital share production function')
    delta       ${\delta}$ (long_name='depreciation rate')
    n_initial   ${n_0}$ (long_name='initial population growth rate')
    g_initial   ${g_0}$ (long_name='initial technology growth rate')
    ;

//****************************************************************************
//set parameters
//****************************************************************************
s=0.2;
alpha=0.3;
delta=0.1;
n_initial=0.01;
g_initial=0.02;
        
//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************
model;
[name='Law of motion capital']
(1+n+g+n*g)*k(+1)=(1-delta)*k+invest;
[name='resource constraint']
invest+c=y;
[name='behavioral rule savings']
c=(1-s)*y;
[name='production function']
y=k^alpha;
[name='Definition log output']
log_y=log(y);
[name='Definition log consumption']
log_c=log(c);
[name='Definition log investment']
log_invest=log(invest);
[name='Definition capital decided upon today']
log_k=log(k(+1));
[name='Definition capital growth rate between today and tomorrow ']
g_k_intensive=log(k(+1))-log(k); //this is capital growth rate between today and tomorrow that is decided upon based on the growth rate of labor augmenting technology between today and tomorrow
[name='Definition capital per capita growth rate between today and tomorrow']
g_k_per_capita=g_k_intensive+g; 
[name='Definition aggregate capital growth rate between today and tomorrow']
g_k_aggregate=g_k_intensive+g+n; 
end;

//****************************************************************************
//initval-block: set initial condition to steady state value with g=g_initial
//and n=n_initial
//****************************************************************************
initval;
    g=g_initial;
    n=n_initial;
    k=((delta+n+g_initial+n*g_initial)/s)^(1/(alpha-1));
    y=k^alpha;
    c=(1-s)*y;
    invest=y-c;
    log_y=log(y);
    log_c=log(c);
    log_invest=log(invest);
    log_k=log(k);
    g_k_intensive=0;
    g_k_per_capita=g_k_intensive+g; 
    g_k_aggregate=g_k_intensive+g+n;
end;

//****************************************************************************
//endval-block: set 
//      i)  g=0 for all periods after the initial one
//      ii) terminal condition to steady state value with g=0;
//****************************************************************************
endval;
    g=0;
    n=n_initial;
    k=((delta+n+g+n*g)/s)^(1/(alpha-1));
    y=k^alpha;
    c=(1-s)*y;
    invest=y-c;
    log_y=log(y);
    log_c=log(c);
    log_invest=log(invest);
    log_k=log(k);
    g_k_intensive=0;
    g_k_per_capita=g_k_intensive+g; 
    g_k_aggregate=g_k_intensive+g+n;
end;
//****************************************************************************
//steady-command: compute the steady state conditional on the value of the exogenous 
//  states and use them for endval (redundant here, because we provided the steady states
//  analytically) 
//****************************************************************************
steady;

//****************************************************************************
//perfect_foresight_setup: set up the simulation for 100 periods
//you can check the settings in oo_.endo_simul (and oo_.exo_simul if there 
//were exogenous variables)
//****************************************************************************
perfect_foresight_setup(periods=100);

//****************************************************************************
//perfect_foresight_solver: compute the solution
//****************************************************************************
perfect_foresight_solver;

//****************************************************************************
//rplot command: display simulation results
//****************************************************************************
rplot k y c invest;

rplot log_k log_y log_c log_invest;

figure('Name','Change in g')
subplot(4,1,1)
plot(log_k)
axis tight
title('Log k in intensive form')
        
subplot(4,1,2)
plot(g_k_intensive)
axis tight
title('Growth rate of k in intensive form')

subplot(4,1,3)
plot(g_k_per_capita)
axis tight
title('Growth rate of per capita k')

subplot(4,1,4)
plot(g_k_aggregate)
axis tight
title('Growth rate of aggregate k')

%%%%%%%%%%%%%% Now consider change in n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//****************************************************************************
//initval-block: set initial condition to steady state value with g=g_initial
//and n=n_initial        
//****************************************************************************
initval;
    g=g_initial;
    n=n_initial;
    k=((delta+n+g_initial+n*g_initial)/s)^(1/(alpha-1));
    y=k^alpha;
    c=(1-s)*y;
    invest=y-c;
    log_y=log(y);
    log_c=log(c);
    log_invest=log(invest);
    log_k=log(k);
    g_k_intensive=0;
    g_k_per_capita=g_k_intensive+g; 
    g_k_aggregate=g_k_intensive+g+n;
end;

//****************************************************************************
//endval-block: set 
//      i)  n=0 for all periods after the initial one
//      ii) terminal condition to steady state value with n=0;
//****************************************************************************
endval;
    g=g_initial;
    n=0;
    k=((delta+n+g+n*g)/s)^(1/(alpha-1));
    y=k^alpha;
    c=(1-s)*y;
    invest=y-c;
    log_y=log(y);
    log_c=log(c);
    log_invest=log(invest);
    log_k=log(k);
    g_k_intensive=0;
    g_k_per_capita=g_k_intensive+g; 
    g_k_aggregate=g_k_intensive+g+n;
end;
//****************************************************************************
//steady-command: compute the steady state conditional on the value of the exogenous 
//  states and use them for endval (redundant here, because we provided the steady states
//  analytically) 
//****************************************************************************
steady;

//****************************************************************************
//perfect_foresight_setup: set up the simulation for 100 periods
//you can check the settings in oo_.endo_simul (and oo_.exo_simul if there 
//were exogenous variables)
//****************************************************************************
perfect_foresight_setup(periods=100);

//****************************************************************************
//perfect_foresight_solver: compute the solution
//****************************************************************************
perfect_foresight_solver;

//****************************************************************************
//rplot command: display simulation results
//****************************************************************************
rplot k y c invest;

rplot log_k log_y log_c log_invest;

figure('Name','Change in n')
subplot(4,1,1)
plot(log_k)
axis tight
title('Log k in intensive form')
        
subplot(4,1,2)
plot(g_k_intensive)
axis tight
title('Growth rate of k in intensive form')

subplot(4,1,3)
plot(g_k_per_capita)
axis tight
title('Growth rate of per capita k')

subplot(4,1,4)
plot(g_k_aggregate)
axis tight
title('Growth rate of aggregate k')
