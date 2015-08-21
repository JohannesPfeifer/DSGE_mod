/*
 * This file replicates the model studied in:
 * Caldara, Dario/Fernandez-Villaverde, Jesus/Rubio-Ramirez, Juan F./Yao, Wen (2012): 
 * "Computing DSGE Models with Recursive Preferences and Stochastic Volatility"
 * Review of Economic Dynamics, 15, pp. 188-206.
 * 
 * It provides a full replication of the main results of the original paper  
 * for second and third order perturbation.
 *
 * This mod-file shows how to use auxiliary variables to deal with recursive preferences
 * and expected returns.
 * 
 * Notes:
 * - The model is written in Dynare's end of period stock notation.
 * - The following things in the original paper are incorrect:
 *      - p.194: the second augmented equilibrium condition is wrong. There is a labor term
 *          (((1-l(+1))/(1-l))^(1-nu))^((1-gamma)/theta) missing
 *      - p. 197: Table 1: the calibrated value for nu is wrong; labor is calibrated to 1/3 and nu chosen
 *          accordingly; the actual value is nu = 0.36218431417051 according to the published 
 *          replication Fortran codes
 *      - p. 197: Table 1: the calibrated value for sigma_bar is wrong; they set sigma_bar=log(0.007) and not
 *          log(sigma_bar)=0.007
 *      - p. 198: Figure 1: the y-label of the right-upper corner should be L(..) not C(..)
 *      - p. 199: Figure 2: the bottom row of graphs seems to be wrong; C has a mean value of about
 *          0.73; it cannot be 0.945 for z=0 and k at the steady state 
 *          (it's almost outside of the density of Figure 4)
 *          Similarly, for labor displayed in the right bottom graph at second order 
 *          should be closer to the mean value for sigma=0.021 shown in Figure 4, which is around 0.33
 *      - p. 201: Table 4: the caption should say "extreme calibration" instead of "benchmark calibration"
 *
 * Requirements: 
 *      - statistics toolbox (uses ksdensity for plotting density estimators)
 *      - function plot_policy_fun available from https://github.com/JohannesPfeifer/DSGE_mod      
 *
 * This implementation was written by Johannes Pfeifer. If you spot any mistakes,
 * email me at jpfeifer@gmx.de
 * Please note that the following copyright notice only applies to this
 * implementation of the model.
 */

/*
 * Copyright (C) 2014-15 Johannes Pfeifer
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

@#define extreme_calibration =1

var V $V$
    y $y$
    c $c$
    k $k$
    invest $i$
    l $l$
    z $z$
    s $s$
    E_t_SDF_plus_1 ${E_t(SDF_{t+1}}$
    sigma $\sigma$
    E_t_R_k ${E_t(R^k_{t+1})}$
    R_f ${R^f}$;

varexo e $\varepsilon$
    omega $\omega$;

parameters beta $\beta$
    gamma $\gamma$
    delta $\delta$
    nu $\nu$
    psi $\psi$
    lambda $\lambda$
    zeta $\zeta$
    rho $\rho$
    sigma_bar ${\bar \sigma}$
    eta $\eta$;

beta = 0.991;
nu = 0.362184314170512;  % typo in paper; fixed nu so that l=1/3 (according to their Fortran codes)
zeta = 0.3;
delta = 0.0196;
lambda = 0.95;

@#if extreme_calibration
    psi = 0.5;
    gamma = 40;
    sigma_bar = log(0.021); % typo in paper; not log(sigma)=0.007
    eta=0.1;
@# else
    psi = 0.5;
    gamma = 5;
    sigma_bar = log(0.007); % typo in paper; not log(sigma)=exp0.007
    eta=0.06;
@# endif

rho=0.9;


model;
#theta = (1 - gamma)/(1 - (1/psi));
// Define Value function
V = ((1-beta)*((c^nu)*((1-l)^(1-nu)))^((1-gamma)/theta) + beta*s^(1/theta))^(theta/(1-gamma));

// Define an auxiliary variable s that captures E_t[V(+1)^sigma]
s =V(+1)^(1-gamma);

// Euler equation: was wrong in paper as labor term was missing
1 = beta*(((1-l(+1))/(1-l))^(1-nu)*(c(+1)/c)^nu)^((1-gamma)/theta)*c/c(+1)
    *((V(+1)^(1-gamma))/s)^(1-(1/theta))*(zeta*exp(z(+1))*k^(zeta-1)*l(+1)^(1-zeta) + 1 - delta); 

//define net return to capital
E_t_R_k=zeta*exp(z(+1))*k^(zeta-1)*l(+1)^(1-zeta) - delta;

//define expected value of stochastic discount factor
E_t_SDF_plus_1=beta*(((1-l(+1))/(1-l))^(1-nu)*(c(+1)/c)^nu)^((1-gamma)/theta)*c/c(+1)
    *((V(+1)^(1-gamma))/s)^(1-(1/theta));

//define net risk-free rate
R_f=(1/E_t_SDF_plus_1-1);
    
// Labor supply FOC
((1-nu)/nu)*(c/(1-l)) = (1-zeta)*exp(z)*(k(-1)^(zeta))*(l^(-zeta)); 

//Budget constraint
c +invest = exp(z)*(k(-1)^(zeta))*(l^(1-zeta));

// Law of motion of capital
k = (1-delta)*k(-1) + invest;

// Technology shock
z = lambda*z(-1) + exp(sigma)*e; 

// Output definition
y = exp(z)*(k(-1)^(zeta))*(l^(1-zeta));

// LOM Volatility
sigma=(1-rho)*sigma_bar+rho*sigma(-1)+eta*omega;
end;


steady_state_model;
%%steady state according to appendix, but not actually used; use this for
%%changing parameter that affect steady state, while keeping nu constant
%Omega=(1/zeta*(1/beta-(1-delta)))^(1/(zeta-1));
%Phi=nu/(1-nu)*(1-zeta)*Omega^zeta;
%l=Phi/(Omega^zeta-delta*Omega+Phi);
%k=Phi*Omega/(Omega^zeta-delta*Omega+Phi);
%c = k^zeta*l^(1-zeta)-delta*k;

%%steady state actually used; sets labor to 1/3 and adjusts nu accordingly
l = 1/3;
k = ((1-beta*(1-delta))/(zeta*beta))^(1/(zeta-1))*l;
c = k^zeta*l^(1-zeta)-delta*k;
nu = c/((1-zeta)*k^zeta*l^(-zeta)*(1-l)+c);
invest = k^zeta*l^(1 - zeta) - c;
V=c^nu*(1-l)^(1-nu);
s = V^(1-gamma);
z = 0;
y=k^zeta*l^(1-zeta);
E_t_SDF_plus_1=beta;
sigma=sigma_bar;
E_t_R_k=zeta*k^(zeta-1)*l^(1-zeta) - delta;
R_f=1/E_t_SDF_plus_1-1;
end;

steady;
check;

shocks;
var e; stderr 1;
var omega; stderr 1;
end;
write_latex_dynamic_model;

%% get second order decision rules
stoch_simul(order=2,periods=100000,drop=1000,irf=0) c l k y E_t_R_k R_f;

%%plot densities
verbatim;
hh=figure;
subplot(3,2,1)
[f,xi]=ksdensity(c);
plot_handle1=plot(xi,f,'r--');
axis tight
@#if extreme_calibration
xlim([0.5 1])
ylim([0 7])
@# else
xlim([0.64 0.8])
ylim([0 20])
@# endif

subplot(3,2,2)
[f,xi]=ksdensity(l);
plot(xi,f,'r--')
axis tight
@#if extreme_calibration
xlim([0.24 0.4])
ylim([0 50])
@# else
xlim([0.31 0.35])
ylim([0 120])
@# endif

subplot(3,2,3)
[f,xi]=ksdensity(k);
plot(xi,f,'r--')
axis tight
@#if extreme_calibration
xlim([6 16])
ylim([0 0.35])
@# else
xlim([8 11])
ylim([0 1.4])
@# endif

subplot(3,2,4)
[f,xi]=ksdensity(y);
plot(xi,f,'r--')
axis tight
@#if extreme_calibration
xlim([0.4 1.8])
ylim([0 5])
@# else
xlim([0.75 1.15])
ylim([0 12])
@# endif

subplot(3,2,5)
[f,xi]=ksdensity(R_f);
plot(xi,f,'r--')
axis tight
@#if extreme_calibration
xlim([-5 20]*10^(-3))
ylim([0 200])
@# else
xlim([5 12]*10^(-3))
ylim([0 600])
@# endif

subplot(3,2,6)
[f,xi]=ksdensity(E_t_R_k);
plot(xi,f,'r--')
axis tight
@#if extreme_calibration
xlim([-5 20]*10^(-3))
ylim([0 200])
@# else
xlim([5 13]*10^(-3))
ylim([0 600])
@# endif
end;


@#if extreme_calibration
k_range=[3:0.1:32];
sigma_range=log([0.01:0.001:0.045]);
@# else
k_range=[5.72:0.01:13.36];
sigma_range=log([4.5:0.01:10.5]*10^(-3));
@# endif


[k_vec_prun_2,c_vec_k_prun_2,policyfun]=plot_policy_fun('k',k_range,'c',oo_.dr.ys,0);
[sigma_vec_prun_2,c_vec_sigma_prun_2,policyfun]=plot_policy_fun('sigma',sigma_range,'c',oo_.dr.ys,0);
[k_vec_prun_2,l_vec_k_prun_2,policyfun]=plot_policy_fun('k',k_range,'l',oo_.dr.ys,0);
[sigma_vec_prun_2,l_vec_sigma_prun_2,policyfun]=plot_policy_fun('sigma',sigma_range,'l',oo_.dr.ys,0);

%% get second order decision rules
stoch_simul(order=3,periods=100000,drop=1000,irf=0);

%%plot densities
verbatim;
figure(hh)
subplot(3,2,1)
hold on
[f,xi]=ksdensity(c);
plot_handle2=plot(xi,f,'b-');
legend([plot_handle1,plot_handle2],'Second Order','Third Order')

subplot(3,2,2)
hold on
[f,xi]=ksdensity(l);
plot(xi,f,'b-')

subplot(3,2,3)
hold on
[f,xi]=ksdensity(k);
plot(xi,f,'b-')

subplot(3,2,4)
hold on
[f,xi]=ksdensity(y);
plot(xi,f,'b-')

subplot(3,2,5)
hold on
[f,xi]=ksdensity(R_f);
plot(xi,f,'b-')

subplot(3,2,6)
hold on
[f,xi]=ksdensity(E_t_R_k);
plot(xi,f,'b-')

end;

%%plot policy functions
[k_vec_prun_3,c_vec_k_prun_3,policyfun]=plot_policy_fun('k',k_range,'c',oo_.dr.ys,0);
[sigma_vec_prun_3,c_vec_sigma_prun_3,policyfun]=plot_policy_fun('sigma',sigma_range,'c',oo_.dr.ys,0);
[k_vec_prun_3,l_vec_k_prun_3,policyfun]=plot_policy_fun('k',k_range,'l',oo_.dr.ys,0);
[sigma_vec_prun_3,l_vec_sigma_prun_3,policyfun]=plot_policy_fun('sigma',sigma_range,'l',oo_.dr.ys,0);

verbatim;
figure
subplot(2,2,1)
plot(k_vec_prun_3,c_vec_k_prun_2,'r--',k_vec_prun_3,c_vec_k_prun_3,'b-')
axis tight
xlabel('K')
ylabel('C(K,z,\sigma)')
legend('Second Order','Third Order')
subplot(2,2,2)
plot(k_vec_prun_3,l_vec_k_prun_2,'r--',k_vec_prun_3,l_vec_k_prun_3,'b-')
axis tight
xlabel('K')
ylabel('L(K,z,\sigma)')
subplot(2,2,3)
plot(exp(sigma_vec_prun_3),c_vec_sigma_prun_2,'r--',exp(sigma_vec_prun_3),c_vec_sigma_prun_3,'b-')
axis tight
temp=get(gca,'YTick');
set(gca,'YTickLabel',num2str(temp','%8.6f'))
xlabel('\sigma')
ylabel('C(K,z,\sigma)')
subplot(2,2,4)
plot(exp(sigma_vec_prun_3),l_vec_sigma_prun_2,'r--',exp(sigma_vec_prun_3),l_vec_sigma_prun_3,'b-')
axis tight
temp=get(gca,'YTick');
set(gca,'YTickLabel',num2str(temp','%8.6f'))
xlabel('\sigma')
ylabel('L(K,z,\sigma)')
end;

