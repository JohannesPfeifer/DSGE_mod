/*
 * This file replicates the open economy model studied in:
 * George McCandless (2008): The ABCs of RBCs - An Introduction to Dynamic 
 *   Macroeconomic Models, Harvard University Press, Chapter 13
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
    pstar       ${P^*}$     (long_name='foreign price level')
    g           $g$         (long_name='growth rate of money stock')
    lambda      $\lambda$   (long_name='TFP')
    b           $B$         (long_name='foreign bonds')
    rf          ${r^f}$     (long_name='foreign interest rate')
    e           $e$         (long_name='exchange rate')
    x           $X$         (long_name='net exports')
    ;

varexo eps_lambda       ${\varepsilon^\lambda}$    (long_name='TFP shock')
       eps_g            ${\varepsilon^g}$          (long_name='Money growth shock')
       eps_pstar        ${\varepsilon^*}$          (long_name='Foreign price level shock')
       ;
       
parameters beta         ${\beta}$    (long_name='discount factor')
           delta        ${\delta}$   (long_name='depreciation rate')
           theta        ${\theta}$   (long_name='capital share production')
           kappa        ${\kappa}$   (long_name='capital adjustment cost')
           a            ${a}$        (long_name='risk premium')
           B            ${B}$        (long_name='composite labor disutility parameter')
           gamma_lambda ${\gamma_\lambda}$  (long_name='autocorrelation TFP')
           gamma_g      ${\gamma_g}$        (long_name='autocorrelation money growth')
           gamma_pstar  ${\gamma_{P^*}}$    (long_name='autocorrelation foreign price')
           pistar       ${\pi^*}$           (long_name='foreign inflation')
           rstar        ${\r^*}$            (long_name='foreign interest rate')
           sigma_lambda ${\sigma_\lambda}$  (long_name='standard deviation TFP shock')
           sigma_g      ${\sigma_g}$        (long_name='standard deviation money shock')
           sigma_pstar  ${\sigma_{P^*}}$    (long_name='standard deviation foreign price shock')
           ;
           

kappa = 0.5;
beta = 0.99;
delta = 0.025;
theta = 0.36;
rstar = 0.03;
a = 0.01;
B = -2.58;
gamma_lambda = 0.95;
gamma_g = 0.95;
gamma_pstar = 0.95;
sigma_lambda = 0.01;
sigma_g = 0.01;
sigma_pstar = 0.01;

model;

[name='Euler equation']
0 = e/(p(+1)*c(+1)) - beta*e(+1)*(1+rf)/(p(+2)*c(+2));
[name='FOC capital']
0 = p/(p(+1)*c(+1))*(1+kappa*(k-k(-1))) - 
    beta*p(+1)/(p(+2)*c(+2))*(r(+1)+(1-delta)+kappa*(k(+1)-k));
[name='FOC hours worked']
0 = B/w + beta*p/(p(+1)*c(+1));
[name='cash in advance constraint']
0 = p*c-m;
[name='Budget constraint']
0 = m/p + e*b/p + k + kappa/2*( k-k(-1) )^2 - 
    w*h - r*k(-1) - (1-delta)*k(-1) - e*(1+rf(-1))*b(-1)/p;
[name='Firm FOC labor']
0 = w - (1-theta)*lambda*k(-1)^theta*h^(-theta);
[name='Firm FOC capital']
0 = r - theta*lambda*k(-1)^(theta-1)*h^(1-theta);
[name='Evolution foreign bonds']
0 = b - (1+rf(-1))*b(-1)-pstar*x;
[name='foreign interest rate']
0 = rf - rstar + a*b/pstar;
[name='Definition exchange rate']
0 = e - p/pstar;
[name='Definition money growth']
0 = m-g*m(-1);
[name='LOM TFP']
lambda = 1-gamma_lambda + gamma_lambda*lambda(-1) + sigma_lambda*eps_lambda;
[name='LOM money']
g = (1-gamma_g)*1 + gamma_g*g(-1) + sigma_g*eps_g;
[name='LOM foreign price']
pstar = (1-gamma_pstar)*1 + gamma_pstar*pstar(-1) + sigma_pstar*eps_pstar;
end;

 
steady_state_model;
pistar = 1;
r = 1/beta-(1-delta);
rf = 1/beta-1;
b = (rstar+1-1/beta)/a;
x = ( (1-beta)^2-(1-beta)*beta*rstar )/(a*beta^2);
w = (1-theta)*(theta/r)^(theta/(1-theta));
c = beta*w/(-B*pistar);
m_pss = c; //m over p 
k = theta*(m_pss-rf*b)/(r-theta*delta);
h =r*(1-theta)/(w*theta)*k;
lambda = 1;
g = 1;
pstar=1;
p = 1;
m = m_pss*p;
e = 1;
end;


shocks;
var eps_lambda; stderr 1;
var eps_g; stderr 1;
var eps_pstar; stderr 1;
end;

resid;

steady;

stoch_simul(order =1, irf=100,periods = 0) k c w b m p e rf r;