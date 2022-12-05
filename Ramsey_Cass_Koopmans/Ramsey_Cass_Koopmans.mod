/*
 * This file uses Dynare's perfect foresight solver to study the transition 
 * behavior of simple non-stationary Ramsey-Cass-Koopmans economy with Cobb-Douglas production function
 * to its balanced growth path (BGP). It starts the simulation with a capital stock 
 * corresponding to 90% of its BGP value.
 * 
 * Notes:
 *  - The Ramsey-Cass-Koopmans model is solved here in aggregate, i.e. non-detrended form. Because
 *      in aggregate form there is no stationary steady state (only a BGP), one cannot
 *      use the steady-command after initval and endval to compute a conditional steady state.
 *  - The initial and terminal conditions are computed by using the analytical steady state
 *      in intensive form and then transform these values back to the aggregate BGP values
 *
 * This implementation was written by Johannes Pfeifer. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2014-2022 Johannes Pfeifer
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

@#define simulation_periods=30

//****************************************************************************
//Define variables
//****************************************************************************
var C       ${C}$ (long_name='consumption')
    K       ${K}$ (long_name='capital')
    Y       ${Y}$ (long_name='output')
    invest  ${I}$ (long_name='investment')
    log_C       ${\log C}$ (long_name='log consumption')
    log_K       ${\log K}$ (long_name='log capital')
    log_Y       ${\log Y}$ (long_name='log output')
    log_invest  ${\log I}$ (long_name='log investment')
    g_K_aggregate       ${\Delta K}$ (long_name='Growth rate aggregate capital')
    g_K_per_capita      ${\Delta \tilde k}$ (long_name='Growth rate capital per capita')
    g_K_intensive       ${\Delta k}$ (long_name='Growth rate capital in intensive form')
    g_Y_aggregate       ${\Delta Y}$ (long_name='Growth rate aggregate output')
    g_Y_per_capita      ${\Delta \tilde Y}$ (long_name='Growth rate output per capita')
    g_Y_intensive       ${\Delta Y}$ (long_name='Growth rate output in intensive form')
        ;

//****************************************************************************
//Define exogenous variables (Labor augmenting technology is growing in our experiment)
//****************************************************************************
varexo  A    ${A}$ (long_name='Labor augmenting technology')
        L    ${L}$ (long_name='Aggregate labor')
    ;

//****************************************************************************
//Define parameters
//****************************************************************************
parameters alpha        ${\alpha}$      (long_name='capital share production function')
    beta                ${\beta}$       (long_name='discount factor')
    delta               ${\delta}$      (long_name='depreciation rate')
    n                   ${n}$           (long_name='population growth rate')
    g                   ${g}$           (long_name='technology growth rate')
    ;

//****************************************************************************
//set parameters
//****************************************************************************

alpha=0.3;
delta=0.1;
beta=0.99;
n=0.01;
g=0.02;

//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************
model;
    [name='Law of motion capital']
    K=(1-delta)*K(-1)+invest;
    [name='resource constraint']
    invest+C=Y;
    [name='behavioral rule savings']
    1/C=beta*1/C(+1)*(alpha*Y(+1)/K+(1-delta));
    [name='production function']
    Y=K(-1)^alpha*(A*L)^(1-alpha);
    [name='Definition log output']
    log_Y=log(Y);
    [name='Definition log consumption']
    log_C=log(C);
    [name='Definition log investment']
    log_invest=log(invest);
    [name='Definition capital decided upon today']
    log_K=log(K);
    [name='Definition aggregate capital growth rate']
    g_K_aggregate=(K-K(-1))/K(-1);
    [name='Definition capital per capita growth rate']
    g_K_per_capita=(K/L-K(-1)/L(-1))/(K(-1)/L(-1)); 
    [name='Definition capital growth rate']
    g_K_intensive=(K/(A*L)-K(-1)/(A(-1)*L(-1)))/(K(-1)/(A(-1)*L(-1))); 
    [name='Definition aggregate output growth rate']
    g_Y_aggregate=(Y-Y(-1))/Y(-1);
    [name='Definition output per capita growth rate']
    g_Y_per_capita=(Y/L-Y(-1)/L(-1))/(Y(-1)/L(-1)); 
    [name='Definition output growth rate']
    g_Y_intensive=(Y/(A*L)-Y(-1)/(A(-1)*L(-1)))/(Y(-1)/(A(-1)*L(-1))); 
end;

//****************************************************************************
//initval-block:  set initial condition to capital of 90% of steady state value 
//and A_0=1 and L_0=1
//****************************************************************************
initval;
    A=1*(1+g)^0; %A_0
    L=1*(1+n)^0; %L_0
    %compute predetermined capital stock, which is decided at time 0, but used for production 
    %at time 1, i.e. K_0 in our notation, which contains A_0 and L_0;
    K=0.9*(A*L)*(1+n)*(1+g)*((1/beta*(1+n)*(1+g)-(1-delta))/(alpha))^(1/(alpha-1)); %steady state capital in intensive form multiplied by A*L
    %Output at time 0 would have been produced with capital determined at
    %time -1, i.e. K(-1) in Dynare's notation; 
    Y=(K/(1+n+g+n*g))^alpha*(A*L)^(1-alpha); %compute Y based on A_0, L_0, and K_(-1)
    invest=(1-(1-delta)/(1+n+g+n*g))*K;
    C=Y-invest;
    log_Y=log(Y);
    log_C=log(C);
    log_invest=log(invest);
    log_K=log(K);
    g_K_aggregate=n+g+n*g;
    g_K_per_capita=g; 
    g_K_intensive=0; 
    g_Y_aggregate=n+g+n*g;
    g_Y_per_capita=g; 
    g_Y_intensive=0; 
end;
check;
//****************************************************************************
//shocks-block: define path of exogenous variables
//****************************************************************************

shock_vals_A=cumprod((1+g)*ones(@{simulation_periods},1))
shock_vals_L=cumprod((1+n)*ones(@{simulation_periods},1))

shocks;
    var A;
    periods 1:@{simulation_periods};
    values (shock_vals_A);
    var L;
    periods 1:@{simulation_periods};
    values (shock_vals_L);
end;

//****************************************************************************
//endval-block: set terminal condition to steady state value
//****************************************************************************
endval;
    A=1*(1+g)^(@{simulation_periods}+1);
    L=1*(1+n)^(@{simulation_periods}+1);
    %compute predetermined capital stock, which is decided at time T+1, i.e. K_(T+1) 
    %in our notation
    K=(A*L)*(1+n)*(1+g)*((1/beta*(1+n)*(1+g)-(1-delta))/(alpha))^(1/(alpha-1)); %steady state capital in intensive form multiplied by A*L
    Y=(K/(1+n+g+n*g))^alpha*(A*L)^(1-alpha); %compute Y_{T+1} based on A_{T+1}, L_{T+1}, and K_{T}
    invest=(1-(1-delta)/(1+n+g+n*g))*K;
    C=Y-invest;
    log_Y=log(Y);
    log_C=log(C);
    log_invest=log(invest);
    log_K=log(K);
    g_K_aggregate=n+g+n*g;
    g_K_per_capita=g; 
    g_K_intensive=0; 
    g_Y_aggregate=n+g+n*g;
    g_Y_per_capita=g; 
    g_Y_intensive=0; 
end;

//****************************************************************************
//perfect_foresight_setup: set up the simulation for 200 periods
//you can check the settings in oo_.endo_simul (and oo_.exo_simul if there 
//were exogenous variables)
//****************************************************************************
perfect_foresight_setup(periods=@{simulation_periods});
//****************************************************************************
//perfect_foresight_solver: compute the solution
//****************************************************************************
perfect_foresight_solver;

//****************************************************************************
//rplot command: display simulation results
//****************************************************************************

rplot log_K log_C log_Y;
rplot g_K_aggregate g_K_per_capita g_K_intensive;
rplot g_Y_aggregate g_Y_per_capita g_Y_intensive;