/*
 * This file implements a simple Susceptible-Infected-Recovered (SIR) model as in 
 * James H. Stock (2020): "Data Gaps and the Policy Response to the Novel Coronavirus". 
 *
 * Notes:
 *  - The model frequency is weekly and simulations are for a year.
 *  - The model itself is purely backward-looking and can efficiently be solved using 
 *      Dynare's perfect foresight solver
 *  - The present model enforces the complementary slackness constraints that the number of
 *      infected cannot be negative and that the number of infected and recovered cannot exceed
 *      the population. This prevents an overshooting in some extreme calibrations.
 *  - The original paper does not clearly specify the process for the infection rate beta; the 
 *      scenarios in this file replicate them qualitatively, but not exactly in terms of quantitative 
 *      results
 *  - A current limitation of Dynare 4.6.1 is that it does not allow for using the LMMCP solver for 
 *      purely backward-looking models. The present mod-file gets around this by specifying a dummy 
 *      forward-looking equation.
 *  - As a consequence, the last period of the plots needs to be ignored. 
 *  - The mod-file starts from a low number of 50 infected. In terms of population percentage, this is
 *      very low and requires the default tolerance to be decreased. Otherwise, 0 infected will be a solution.
 *
 *
 * THIS MOD-FILE REQUIRES DYNARE 4.6 OR HIGHER
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model. 
 */

/*
 * Copyright (C) 2020 Johannes Pfeifer
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

 
@#ifndef Scenario_A
    @#define Scenario_A=0
@#endif

@#ifndef Scenario_B
    @#define Scenario_B=1
@#endif

@#ifndef Scenario_C
    @#define Scenario_C=0
@#endif

var S                   $S$                     (long_name='population share of susceptible')
    I                   $I$                     (long_name='population share of infected')
    R                   $R$                     (long_name='population share of recovered')
    P_I_given_symptoms  ${P(I|Symptomatic)}$    (long_name='probability of infection given symptoms')
    P_symptomatic       ${P(Symptomatic)}$      (long_name='probability of being symptomatic')
    R_0                 ${R_0)}$                (long_name='basic reproduction number')
    Total_infected      ${I+R}$                 (long_name='total number of infected')
    P_I_and_Symptoms    ${P(I and S)}$          (long_name='Share of infected and symptomatic')
    junk                ${junk}$                (long_name='dummy forward-looking variable')
;

varexo beta
;

parameters
N_total     ${N}$               (long_name='Total population')
gamma       ${\gamma}$          (long_name='recovery rate')
pi_a        ${\pi_a}$           (long_name='share of asymptomatic')
s_0         ${s_0}$             (long_name='baseline rate of symptoms')
;

gamma=0.55;
pi_a=0.3;
s_0=0.02;
N_total=329000000;

model;
[mcp='S>0',name='LOM susceptible']
S-S(-1)=-beta*I(-1)*S(-1);
[mcp='R<1',name='LOM recovered']
R-R(-1)=gamma*I(-1);
[mcp='I<1',name='LOM infected']
I-I(-1)=beta*I(-1)*S(-1)-gamma*I(-1);
[name='Definition total number of people ever infected']
P_I_given_symptoms=(1-pi_a)*I/P_symptomatic;
[name='Share of people showing symptoms']
P_symptomatic=P_I_and_Symptoms+s_0*(S+R);
[name='Definition number of people infected and symptomatic']
P_I_and_Symptoms=(1-pi_a)*I;
[name='Definition total number of people ever infected']
Total_infected=I+R;
[name='Definition basic reproduction number']
R_0=beta/gamma;
[name='dummy equation to work around Dynare limitation']
junk=junk(+1);
end;

initval;
I=50/N_total;
R=0/N_total;
S=1-I-R;
beta=2.1;
P_I_and_Symptoms=(1-pi_a)*I;
P_symptomatic=P_I_and_Symptoms+s_0*(S+R);
R_0=beta/gamma;
end;
resid;

endval;
I=0;
R=1;
S=0;
beta=2.1;
P_I_and_Symptoms=(1-pi_a)*I;
P_symptomatic=P_I_and_Symptoms+s_0*(S+R);
R_0=beta/gamma;
Total_infected=I+R;
end;

@#if Scenario_A
shocks;
    var beta;
    periods 1:12, 13:17, 18,      19,     20,     21,     22,    23,       24:52;
    values 2.1,    0.99 , 1.1750, 1.1750, 1.5450, 1.5450, 1.9150 ,1.9150,  2.1 ;
end;
@#endif

@#if Scenario_B
shocks;
    var beta;
    periods 1:12, 13:24, 25:28,   29:32, 33:36 ,  37:40,   41:44, 45:52;
    values  2.1,  0.66 , 0.8800, 0.9900, 1.1000,    1.2100,    1.3200, 2.1 ;
end;
@#endif


@#if Scenario_C
shocks;
    var beta;
    periods 1:12, 13:17,17:52;
    values 2.1, 0.1760 ,0.4;
end;
@#endif

perfect_foresight_setup(periods=52);

% get solution neglecting natural bounds on R,I,S and use that as a starting point
% Attention: the tolerance level needs to be low, because the starting
% value of 50 infected people is in the range of 1e-7
perfect_foresight_solver(tolf=1e-8);

perfect_foresight_solver(lmmcp); 

rplot R_0;
rplot I Total_infected;
rplot I R S;
rplot P_I_and_Symptoms;
