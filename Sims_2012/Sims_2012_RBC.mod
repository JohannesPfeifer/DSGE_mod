/*
 * This file replicates the results for the basic RBC model presented in Eric R. Sims (2012):
 * "New, Non-Invertibility, and Structural VARs", Advances in Econometrics, Volume 28, 81–135
 * 
 * In particular, it performs the ABCD-test of Fernandez-Villaverde, Rubio-Ramirez,Sargent, 
 * and Watson (2007), "ABCs (and Ds) of Understanding VARs", American
 * Economic Review, 97(3), 1021-1026

 * This mod-file strips away all features from the full model not used in the RBC model.
 * In particular, epsilon has been set to Infinity, and investment adjustment costs, habits, and 
 * monetary policy have been dropped from the model.
 *
 * This file require the ABCD_test.m from the FV_et_al_2007-folder to be located in the same folder.
 *
 * This file was written by Johannes Pfeifer.In case you spot mistakes, email me at jpfeifer@gmx.de
 *
 * The model is written in the beginning of period stock notation. To make the model
 * conform with Dynare's end of period stock notation, I use the 
 * predetermined_variables-command.
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright (C) 201315 Johannes Pfeifer
 *
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can receive a copy of the GNU General Public License
 * at <http://www.gnu.org/licenses/>.
 */

var c lambda $\lambda$ w n R y mu_a mu_y invest k z1 z2 z3; //interest pi //drop monetary policy
varexo epsilon u;

predetermined_variables k;

parameters  alpha $\alpha$ 
    beta $\beta$ 
    theta $\theta$
    delta $\delta$
    xi $\xi$
    delta_I ${\Delta_I}$    
    rho $\rho$
    phi_pi $\phi_\pi$
    phi_y $\phi_y$
    pi_star $pi^*$
    delta_y ${\Delta_y}$
    g_a $g_a$
    sigma_eps $\sigma_eps$
    sigma_u $\sigma_u$;
    
    
beta=0.99;
delta=0.02;
pi_star=0.005;
alpha=1/3;
xi=1;
sigma_eps=0.01;
sigma_u=0.005;
theta=0; //set in SS
rho=0.8;
phi_pi=1.5;
phi_y=1.5;
g_a=0.0025;
delta_y=exp(g_a)^(1/(1-alpha));
delta_I=delta_y;

model;
//1. Resource Constraint
c+invest=exp(y);
//2. LOM capital
exp(k(+1))=invest+(1-delta)*exp(k)/mu_y;
//3. Lagrange multiplier
lambda=1/c;
//4. FOC labor
theta*n^xi=lambda*w;
//5. Euler equation Bonds
//lambda=beta/mu_y(+1)*lambda(+1)*(1+interest)/(1+pi(+1));
//6. Euler equation capital
lambda=beta/mu_y(+1)*lambda(+1)*(R(+1)+(1-delta));
//7. Output
exp(y)=(exp(k)/mu_y)^alpha*n^(1-alpha);
//8. FOC labor
w=(1-alpha)*exp(y)/n;
//9. FOC Interest
R=alpha*exp(y)/(exp(k)/mu_y);
//10. Taylor Rule
//interest=rho*interest(-1)+(1-rho)*steady_state(interest)+phi_pi*(pi-pi_star)+(1-rho)*phi_y*(y/y(-1)-1);
mu_a=g_a+epsilon+z1(-1); //mu_a=ln a_t - ln a_{t-1}
mu_y=exp(mu_a)^(1/(1-alpha));
z1=z2(-1);
z2=z3(-1);
z3=u;
end;

shocks;
var epsilon; stderr sigma_eps;
var u; stderr sigma_u;
end;

steady_state_model;
mu_a=g_a;
mu_y=exp(mu_a)^(1/(1-alpha));
n=1/3; //calibrate 
k_bar=((1/(beta/mu_y)-(1-delta))/(alpha))^(1/(alpha-1))*n*mu_y;
w=(1-alpha)*(k_bar/mu_y/n)^alpha;
R=alpha*(k_bar/mu_y/n)^(alpha-1);
invest=(1-(1-delta)/mu_y)*k_bar;
pi=pi_star;
y_bar=(k_bar/mu_y)^alpha*n^(1-alpha);
c=y_bar-invest;
lambda=1/c;
interest=1/(beta/mu_y)*(1+pi)-1;
theta=lambda*w/n^xi;
k=log(k_bar);
y=log(y_bar);
end;
steady;
check;

write_latex_dynamic_model;

varobs mu_a y;

stoch_simul(order=1,irf=20) k y mu_a mu_y;


//************************************* Replicate the theoretical IRFs of Figure 2: Monte Carlo Results in RBC Model
//Now back out IRF of non-stationary model variables by adding trend growth back

log_a_surprise(1,1)=0; 
log_a_antic(1,1)=0; 

log_mu_y_surprise_0=0; //Initialize Level of Technology at t=0;
log_mu_y_surprise(1,1)=log_mu_y_surprise_0; 
log_y_nonstationary_surprise(1,1)=0;

log_mu_y_antic_0=0; //Initialize Level of Technology at t=0;
log_mu_y_antic(1,1)=log_mu_y_antic_0; 
log_y_nonstationary_antic(1,1)=0;

// reaccumulate the non-stationary level series
for ii=2:options_.irf+1
    log_a_surprise(ii,1)=log_a_surprise(ii-1,1)+mu_a_epsilon(ii-1,1);
    log_a_antic(ii,1)=log_a_antic(ii-1,1)+mu_a_u(ii-1,1);
    log_mu_y_surprise(ii,1)=mu_y_epsilon(ii-1,1)+log_mu_y_surprise(ii-1,1);
    log_y_nonstationary_surprise(ii,1)=y_epsilon(ii-1,1)+log_mu_y_surprise(ii,1);
    log_mu_y_antic(ii,1)=mu_y_u(ii-1,1)+log_mu_y_antic(ii-1,1);
    log_y_nonstationary_antic(ii,1)=y_u(ii-1,1)+log_mu_y_antic(ii,1);
end

//Make the plot
figure('Name','IRF of non-detrended variables to a productivity shock');
subplot(4,1,1)
plot(1:options_.irf+1,log_a_surprise)
title('Surprise TFP')
subplot(4,1,2)
plot(1:options_.irf+1,log_y_nonstationary_surprise)
title('Output Surprise')
subplot(4,1,3)
plot(1:options_.irf+1,log_a_antic)
title('Antic TFP')
subplot(4,1,4)
plot(1:options_.irf+1,log_y_nonstationary_antic)
title('Antic Output')


//*************************************perform ABCD-test on model, replicates Appendix A.2

[result,eigenvalue_modulo,A,B,C,D]=ABCD_test(M_,options_,oo_)

