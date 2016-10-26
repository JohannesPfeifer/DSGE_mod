/*
 * This file implements the deterministic optimal policy exercise in Figure 7.1 of Michael Woodford (2003): "Interest and prices",
 * Princeton University Press, page 476. The same figure is reproduced as Figure 2 in Michael Woodford (2010): "Optimal Monetary Stabilization
 *  Policy", Chapter 14 of the Handbook of Monetary Economics, Volume 3B, Elsevier
 *
 * The mod-file demonstrates how to use the ramsey_model-command to let Dynare generate the stationary Ramsey first order conditions
 * in order to conduct Ramsey policy from a timeless perspective (set Ramsey_policy_timeless=1). If Ramsey_policy_t0_optimal
 * it uses the t_0 optimal Ramsey policy where the initial multiplier is set to 0.
 * 
 * THIS MOD-FILE REQUIRES DYNARE 4.5 (I.E. THE CURRENT UNSTABLE VERSION)
 *
 * Notes:
 *  - In contrast to what is stated in the article, the parameterization used for Figure 2 is not the one given in Woodford (2003).
 *      While Figure 7.1 is identical to Figure 2, the parameterization is the one used in Woodford (1999), i.e. lambda=0.003 instead of 0.048
 *  - The inflation rate in the Figure is the annualized inflation rate
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



var x       ${x}$ (long_name='output gap')
    pi      ${\pi}$ (long_name='quarterly net inflation rate')
    pi_annual ${\pi^{ann}}$ (long_name='annualized inflation in percent')
            ;

@#define Ramsey_policy_timeless=1
@#define Ramsey_policy_t0_optimal=0
@#define discretionary_policy=0

parameters lambda   ${\lambda}$ (long_name='output gap')
        x_star      ${x^*}$ (long_name='output gap')
        betta       ${\beta}$ (long_name='output gap')
        kappa       ${\kappa}$ (long_name='output gap')
        theta       ${\theta}$ (long_name='output gap')
        ;

x_star=0.2;     %according to footnote 10, p. 475
betta=0.99;     %Table 6.1, p. 431
kappa=0.024;    %Table 6.1, p. 431  
theta=7.88;     %taken from Woodford (1999)


@#if Ramsey_policy_timeless
    model(linear);
        [name='New Keynesian Phillips Curve, Eq. 8.1.1']
        pi-kappa*x-betta*pi(+1)=0;
        [name='Definition anualized inflation']
        pi_annual=400*pi;
    end;
    steady_state_model;
        lambda=kappa/theta; %taken from Woodford (1999) instead of lambda_x from Table 6.1, p. 431 of Woodford (2003)
    end;

    planner_objective(0.5*(pi^2+lambda*(x-x_star)^2));
    ramsey_model(planner_discount=0.99);
    steady;
    simul(periods=50);
@#else
    @#if  Ramsey_policy_t0_optimal
        var varphi;
        model(linear);
            [name='New Keynesian Phillips Curve, Eq. 7.1.1']
            pi-kappa*x-betta*pi(+1)=0;
            [name='FOC 1, Eq. 7.1.7']
            pi+varphi-varphi(-1)=0;
            [name='FOC 2, Eq. 7.1.8']
            lambda*(x-x_star)-kappa*varphi=0;
            [name='Definition anualized inflation']
            pi_annual=400*pi;
        end;

        steady_state_model;
            lambda=kappa/theta; %taken from Woodford (1999) instead of lambda_x from Table 6.1, p. 431 of Woodford (2003)
            x=0;
            pi=0;
            varphi=-lambda*x_star/kappa;
        end;

        histval;
            varphi(0)=0;
        end;

        steady;
        simul(periods=50);        
    @#else    
        model(linear);
            [name='FOC, equation 7.1.4']
            pi+lambda/kappa*(x-x_star)=0;
            [name='equation 7.1.3 with pi=pi^e imposed']
            pi=kappa*x+betta*pi;
            [name='Definition anualized inflation']
            pi_annual=400*pi;
        end;

        steady_state_model;
            lambda=kappa/theta; %taken from Woodford (1999) instead of lambda_x from Table 6.1, p. 431 of Woodford (2003)
            pi=kappa*lambda/((1-betta)*lambda+kappa^2)*x_star;
            x=(1-betta)*pi/kappa;
            pi_annual=400*pi;
        end;
        steady;
        simul(periods=50);
    @#endif
@#endif
rplot pi_annual;
