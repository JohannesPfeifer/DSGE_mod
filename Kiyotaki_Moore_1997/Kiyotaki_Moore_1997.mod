/*
 * This file simulates a version of the model studied in Section 2 of:
 * Kiyotaki/Moore (1997): "Credit Cycles", 
 * Journal of Political Economy, 105(2), 211-248
 * 
 * The implementation is based on Eric Sims lecture notes for his 2020 course
 * "Advanced Macroeconomics", available at https://sites.nd.edu/esims/files/2023/05/kiyotaki_moore_ers_final.pdf 
 *
 * This implementation was written by Johannes Pfeifer based on the codes of Eric Sims. 
 * It uses the production function G=(z+k')^alpha.
 *
 * If you spot mistakes, email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2025 Johannes Pfeifer
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


var x   ${x}$           (long_name='consumption farmer')
    xp  ${x^\prime}$    (long_name='consumption gatherer')
    b   $b$             (long_name='bonds')
    k   ${k}$           (long_name='Land held by farmer')
    kp  ${k^\prime}$    (long_name='Land held by gatherer')
    q   $q$             (long_name='fruit price')
    mu  $\mu$           (long_name='multiplier enforceability constraint')  
    phi $\phi$          (long_name='multiplier consumption constraint farmer')  
    C   $C$             (long_name='aggregate consumption')
    Y   $Y$             (long_name='aggregate output')
    ;

varexo ed $\varepsilon$;

parameters alpha $\alpha$           (long_name='exponent production function gatherers')
    m           $m$                 (long_name='population size gatherers')
    K_bar       $\bar K$            (long_name='land supply')
    beta        $\beta$             (long_name='discount factor farmer')
    betap       ${\beta^\prime}$    (long_name='discount factor farmer')
    a           $a$                 (long_name='tradability share')
    c           $c$                 (long_name='non-tradability share') 
    z           $z$                 (long_name='constant production function gatherers')
    ;

alpha = 1/3;
m = 0.5;
K_bar = 1;
betap = 0.99;
beta = 0.98;
a = 0.7;
c = 0.3;
z = 0.01;

model;
[name='Euler equation bonds farmer']
1 + phi = (beta*(1+phi(+1)) + mu)/betap;

[name='Euler equation capital farmer']
q*(1+phi) + beta*c*phi(+1) = beta*(1+phi(+1))*((1+ed(+1))*(a + c) + q(+1)) + mu*q(+1);

[name='Budget constraint']
q*(k - k(-1)) + b(-1)/betap + x = (1+ed)*(a+c)*k(-1) + b;

[name='Borrowing constraint']
b = betap*q(+1)*k;

[name='Euler equation gatherer']
q = betap*((1+ed(+1))*alpha*(z + kp)^(alpha-1) + q(+1));

[name='Resource constraint']
x + m*xp = (1+ed)*(a+c)*k(-1) + m*(1+ed)*(z + kp(-1))^(alpha);

[name='Capital market-clearing']
k + m*kp = K_bar;

[name='non-tradeable constraint']
x = c*k(-1);

[name='Aggregate consumption']
C = x + m*xp;

[name='Aggregate output']
Y = C;
end;

steady_state_model;
q = a/(1-betap);
kp = (betap*alpha/a)^(1/(1-alpha)) - z;
k = K_bar - m*kp;
b = betap*q*k;
xp = (1/m)*(a*k + m*(z + kp)^(alpha));
phi = (a*(beta-1) + beta*c)/(a*(1-beta));
mu = (betap-beta)*beta*c/(a*(1-beta));
x = c*k;
C = x + m*xp;
Y = C;
end;

shocks;
var ed = 0.0011^2;
end;

stoch_simul(order=1,irf=12,ar=0,TeX) k kp Y q mu;
