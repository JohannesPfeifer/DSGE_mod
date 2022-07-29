/*
 * This file replicates the Money in the Utility Function model studied in:
 * George McCandless (2008): The ABCs of RBCs - An Introduction to Dynamic 
 *   Macroeconomic Models, Harvard University Press, Chapter 9
 * 
 * This implementation was written by Johannes Pfeifer.
 *
 * If you spot mistakes, email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright © 2022 Johannes Pfeifer
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
 * see <https://www.gnu.org/licenses/>.
 */
 
var w           $W$         (long_name='real wage')
    r           $r$         (long_name='real return on capital')
    c           $C$         (long_name='real consumption')
    k           $K$         (long_name='capital stock')
    h           $H$         (long_name='hours worked')
    m           $M$         (long_name='money stock')
    p           $P$         (long_name='price level')
    g           $g$         (long_name='growth rate of money stock')
    lambda      $\lambda$   (long_name='TFP')
    y           $y$         (long_name='real output')
    ;

varexo eps_lambda       ${\varepsilon^\lambda}$    (long_name='TFP shock')
       eps_g            ${\varepsilon^g}$          (long_name='Money growth shock')
       ;
       
parameters beta         ${\beta}$    (long_name='discount factor')
           delta        ${\delta}$   (long_name='depreciation rate')
           theta        ${\theta}$   (long_name='capital share production')
           A            ${A}$        (long_name='labor disutility parameter')
           h_0          ${h_0}$      (long_name='steady state hours worked')
           B            ${B}$        (long_name='composite labor disutility parameter')
           gamma        ${\gamma}$   (long_name='autocorrelation TFP')
           pi           ${\pi}$      (long_name='autocorrelation money growth')
           g_bar        ${\bar g}$   (long_name='steady state growth rate of money')
           D            ${D}$        (long_name='coefficient log balances')
           ;
           
predetermined_variables k;

%Calibration based on Table 9.1
beta = 0.99;
delta = 0.025;
theta = 0.36;
A = 1.72;
h_0 = 0.583;
gamma = 0.95;
pi = 0.48;
g_bar = 1;
D = 0.01;

model;
%composite labor parameter
[name='Budget constraint, (9.1)']
c+k(+1)+m/p = w*h+r*k+(1-delta)*k+m(-1)/p + (g-1)*m(-1)/p;
[name='Euler equation, (9.2)']
1/c = beta*p/(c(+1)*p(+1)) + D*p/m;
[name='FOC hours worked, (9.3)']
1/c = -B/w;
[name='FOC capital, (9.4)']
1/c = beta/c(+1)*(r(+1)+1-delta);
[name='Definition money growth, before (9.5)']
m = g*m(-1);
[name='Production function, below (9.5)']
y = lambda*k^theta*h^(1-theta);
[name='Firm FOC labor, below (9.5)']
w = (1-theta)*lambda*k^theta*h^(-theta);
[name='Firm FOC capital, below (9.5)']
r = theta*lambda*(k/h)^(theta-1);
[name='Law of motion money stock, below (9.5)']
log(g) = (1-pi)*log(g_bar) + pi*log(g(-1))+eps_g;
[name='Law of motion technology shock,  below (9.5)']
log(lambda) = gamma*log(lambda(-1)) + eps_lambda;
end;


steady_state_model;
%follows Section 9.2
B = A*log(1-h_0)/h_0;
r = 1/beta -1 + delta;
w = (1-theta)*(r/theta)^(theta/(theta-1));
c = -w/B;
mp = D*g_bar*c/(g_bar-beta);
k = c/((r*(1-theta)/(w*theta))^(1-theta)-delta);
h = r*(1-theta)/(w*theta)*k;
y = k*(r*(1-theta)/(w*theta))^(1-theta);
g = 1;
lambda = 1;
p = 1; %normalization
m=p*D*g*c/(g-beta);
end;

steady;

shocks;
var eps_g;
stderr 0.01;
end;

stoch_simul(irf=100, order=1) k c w r h m y g p;

shocks(overwrite);
var eps_lambda;
stderr 0.01;
end;

stoch_simul(irf=100, order=1) k c w r h m y g p;
