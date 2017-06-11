/*
 * This file estimates the New Keynesian model of Peter Ireland (2004): Technology shocks in the New Keynesian Model,
 * Review of Economics and Statistics, 86(4), pp. 923-936
 *
 * The mod-file also replicates the variance decompositions and IRFs for the respective samples.
 *
 * Notes: 
 *      - The mod-file requires Dynare 4.5 or higher. Otherwise, you will get error like 
 *          "Check whether your model in truly linear" 
 *      - The original paper is not explicit how the constant in the observation equations are treated during estimation.
 *          It turns out that Ireland individally demeans the observables for each sample so that the observation equations
 *          do not feature a constant. In contrast to what page 928 seems to say, it is not the case that full sample growth 
 *          rates are used for demeaning and only then the sample is cut.
 *      - Ireland restricts omega, alpha_x, and alpha_p to the interval [0,1] by parameter transformation. In this mod-file, the 
 *          constraint is imposed by providing parameter bounds in the estimated_params-block
 *      - Some of the estimated parameters in Ireland are extremely close to the lower bound of 0. This creates problems in the computation
 *          of the Hessian at the mode as two-sided numerical derivatives cannot be computed. Ireland (2004) sidesteps this issue by using 
 *          one-sided derivatives and an individually selected step size. Dynare does not (yet) offer this tailored approach as it is not generic.
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

@#define pre_1980=0
@#define post_1980=1
@#define full_sample=0
 
var a ${a}$ (long_name='preference shock')
    e ${e}$ (long_name='cost push shock')
    z ${az}$ (long_name='TFP shock')
    x ${x}$ (long_name='output gap')
    pihat ${\hat p}$ (long_name='inflation deviation from trend')
    yhat ${\hat y}$ (long_name='output deviations from trend')
    ghat ${\hat g}$ (long_name='output growth')
    rhat ${\hat r}$ (long_name='interest deviations from trend')
    gobs ${g^{obs}}$ (long_name='observed output growth')
    robs ${r^{obs}}$ (long_name='observed interest deviations from trend')
    piobs ${\pi^{obs}}$ (long_name='observed inflation deviations from trend')
    r_annual ${r^{ann}}$ (long_name='annualized interest rate')
    pi_annual ${\pi^{ann}}$ (long_name='annualized inflation rate')
        ; 

varexo eps_a ${\varepsilon_a}$ (long_name='preference innovation')
        eps_e ${\varepsilon_e}$ (long_name='cost push innovation')
        eps_z ${\varepsilon_z}$ (long_name='TFP innovation')
        eps_r ${\varepsilon_r}$ (long_name='monetary policy innovation')
        ;

%% Parameter decleration
parameters beta ${\beta}$ (long_name= 'discount factor')
        alpha_x ${\alpha}$
        alpha_pi ${\alpha_\pi}$
        rho_a ${\rho_a}$
        rho_e ${\rho_e}$
        omega ${\omega}$
        psi ${\psi}$
        rho_pi ${\rho_\pi}$ (long_name= 'inflation feedback')
        rho_g ${\rho_g}$ (long_name= 'output growth feedback')
        rho_x ${\rho_x}$ (long_name= 'output gap feedback')
        ;

%calibrated 
beta = 0.99;    % Ireland's calibration
psi = 0.1;      % Ireland's calibration

% %set parameters to estimated values
@#if full_sample==1 
    omega=0.0617;
    alpha_x = 0.0836;
    alpha_pi = 0.0001;
    rho_pi = 0.3597;
    rho_g = 0.2536;
    rho_x= 0.0347;
    rho_a= 0.9470;
    rho_e= 0.9625;
@#endif

@#if pre_1980==1 
    omega=0.00001;
    alpha_x = 0.2028;
    alpha_pi = 0.00001;
    rho_pi = 0.3053;
    rho_g = 0.2365;
    rho_x= 0.00001;
    rho_a= 0.9910;
    rho_e= 0.5439;
@#endif

@#if post_1980==1 
    omega=0.0581;
    alpha_x = 0.00001;
    alpha_pi = 0.00001;
    rho_pi = 0.3866;
    rho_g = 0.3960;
    rho_x= 0.1654;
    rho_a= 0.9048;
    rho_e= 0.9907;
@#endif

%%   Model Declaration
model(linear);
[tag='temporary preference shock (15)']
a=rho_a*a(-1)+eps_a; 
[tag='temporary cost-push shock (16)']
e=rho_e*e(-1)+eps_e; 
[tag='technology shock (17)']
z=eps_z;
[tag='New Keynesian IS curve (23)']
x=alpha_x*x(-1)+(1-alpha_x)*x(+1)-(rhat-pihat(+1))+(1-omega)*(1-rho_a)*a;
[tag='New Keynesian PC curve (24)']
pihat=beta*(alpha_pi*pihat(-1)+(1-alpha_pi)*pihat(+1))+psi*x-e; 
[tag='output gap (20)']
x=yhat-omega*a;
[tag='growth rate of output (21)']
ghat=yhat-yhat(-1)+z; 
[tag='policy rule (22)']
rhat-rhat(-1)=rho_pi*pihat+rho_g*ghat+rho_x*x+eps_r;
[tag='observation equation g']
gobs=ghat;
[tag='observation equation r']
robs=rhat;
[tag='observation equation pi']
piobs=pihat;
[tag='definition annualized interest rate']
r_annual=4*rhat;
[tag='definition annualized inflation rate']
pi_annual=4*pihat;        
end;


shocks;
@#if full_sample ==1 
    var eps_a; stderr 0.0405;
    var eps_e; stderr 0.0012;
    var eps_z; stderr 0.0109;
    var eps_r; stderr 0.0031;
@#endif
@#if pre_1980 ==1 
    var eps_a; stderr 0.1538;
    var eps_e; stderr 0.0035;
    var eps_z; stderr 0.0104;
    var eps_r; stderr 0.0033;
@#endif
@#if post_1980 ==1 
    var eps_a; stderr 0.0302;
    var eps_e; stderr 0.0002;
    var eps_z; stderr 0.0089;
    var eps_r; stderr 0.0028;
@#endif
end;

% Estimation
estimated_params;
omega; 
alpha_x, ,0,1;
alpha_pi, ,0,1;
rho_pi, ,0,1;
rho_g, ,0,1;
rho_x, ,0,1;
rho_a, ,0,1;
rho_e, ,0,1;
stderr eps_a, ,0,1; 
stderr eps_e, ,0,1;
stderr eps_z, ,0,1;
stderr eps_r, ,0,1;
end;

estimated_params_init(use_calibration);
end;

varobs gobs robs piobs;

% @#if full_sample ==1 
% 	estimation(datafile=data_full_sample,mode_compute=4);
% @#endif
% @#if pre_1980 ==1 
% 	estimation(datafile=data_pre_1980,mode_compute=9);
% @#endif
% @#if post_1980 ==1 
% 	estimation(datafile=data_post_1980,mode_compute=9);
% @#endif
% 
stoch_simul(order=1,conditional_variance_decomposition=[1 4 8 12 20 40], irf=16) ghat pi_annual r_annual x;

figure
subplot(4,4,1)
plot([0:options_.irf],[0 oo_.irfs.ghat_eps_a]*100) 
axis tight
ylabel('Output growth')
title('Preference Shock')

subplot(4,4,5)
plot([0:options_.irf],[0 oo_.irfs.pi_annual_eps_a]*100) 
axis tight
ylabel('Inflation')

subplot(4,4,9)
plot([0:options_.irf],[0 oo_.irfs.r_annual_eps_a]*100) 
axis tight
ylabel('Interest rate')
        
subplot(4,4,13)
plot([0:options_.irf],[0 oo_.irfs.x_eps_a]*100) 
axis tight
ylabel('Output gap')
        
%Cost push shock
subplot(4,4,2)
plot([0:options_.irf],[0 oo_.irfs.ghat_eps_e]*100) 
axis tight
title('Cost push shock')

subplot(4,4,6)
plot([0:options_.irf],[0 oo_.irfs.pi_annual_eps_e]*100) 
axis tight

subplot(4,4,10)
plot([0:options_.irf],[0 oo_.irfs.r_annual_eps_e]*100) 
axis tight
        
subplot(4,4,14)
plot([0:options_.irf],[0 oo_.irfs.x_eps_e]*100) 
axis tight
        
%Technology shock
subplot(4,4,3)
plot([0:options_.irf],[0 oo_.irfs.ghat_eps_z]*100) 
axis tight
title('Technology shock')

subplot(4,4,7)
plot([0:options_.irf],[0 oo_.irfs.pi_annual_eps_z]*100) 
axis tight

subplot(4,4,11)
plot([0:options_.irf],[0 oo_.irfs.r_annual_eps_z]*100) 
axis tight
        
subplot(4,4,15)
plot([0:options_.irf],[0 oo_.irfs.x_eps_z]*100) 
axis tight

%Monetary policy shock
subplot(4,4,4)
plot([0:options_.irf],[0 oo_.irfs.ghat_eps_r]*100) 
axis tight
title('Monetary policy shock')

subplot(4,4,8)
plot([0:options_.irf],[0 oo_.irfs.pi_annual_eps_r]*100) 
axis tight

subplot(4,4,12)
plot([0:options_.irf],[0 oo_.irfs.r_annual_eps_r]*100) 
axis tight
        
subplot(4,4,16)
plot([0:options_.irf],[0 oo_.irfs.x_eps_r]*100) 
axis tight
