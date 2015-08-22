/*
 * This file replicates the neoclassical growth model for Schmitt-Grohé/Uribe (2004): 
 * "Solving dynamic general equilibrium models using a second-order approximation to 
 * the policy function", Journal of Economic Dynamics & Control, 28, pp. 755 – 775
 * 
 * It generates the policy functions in section 5.1. The Dynare output is given by: 
 *
 * % POLICY AND TRANSITION FUNCTIONS
 * %                                    c               k               a
 * % Constant                   -0.969516       -1.552215               0
 * % (correction)               -0.096072        0.241022               0
 * % k(-1)                       0.252523        0.419109               0
 * % epsilon                     0.841743        1.397031        1.000000
 * % k(-1),k(-1)                -0.002559       -0.003501               0
 * % epsilon,epsilon            -0.028433       -0.038901               0
 * % k(-1),epsilon              -0.017060       -0.023341               0
 *
 * Note that the response to the shock epsilon is what SGU call \hat A_t due to no persistence. 
 * Moreover, the second order terms in Dynare are presented including the factor 1/2. 
 * For example the response of consumption to capital^2 
 * is given by 1/2*(-0.0051) in SGU, which equals the -0.002559 provided by Dynare above.
 * 
 * This file was written by Johannes Pfeifer as a Dynare adaptation of the original neoclassical_model.m 
 * provided by SGU in their toolkit. In case you spot mistakes, email me at jpfeifer@gmx.de
 *
 * The model is written in Dynare's end of period stock notation.
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
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

var c k a;
varexo epsilon;

predetermined_variables k;

parameters 	SIG DELTA ALFA BETTA RHO;
BETTA=0.95; %discount rate
DELTA=1; %depreciation rate
ALFA=0.3; %capital share
RHO=0; %persistence of technology shock
SIG=2; %intertemporal elasticity of substitution

model;
0 = exp(c) + exp(k(+1)) - (1-DELTA) * exp(k) - exp(a) * exp(k)^ALFA;
0 = exp(c)^(-SIG) - BETTA * exp(c(+1))^(-SIG) * (exp(a(+1)) * ALFA * exp(k(+1))^(ALFA-1) + 1 - DELTA);
0 = a - RHO * a(-1)-epsilon;
end;

steady_state_model;
k = log(((1/BETTA+DELTA-1)/ALFA)^(1/(ALFA-1))); 
c = log(exp(k)^(ALFA)-DELTA*exp(k)); %steady-state value of consumption 
a = 0; 
end;

shocks;
var epsilon; stderr 1; //from eta=[0 1]'; %Matrix defining driving force
end;

steady;
check;

stoch_simul(order=2);


