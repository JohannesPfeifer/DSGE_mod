/*
 * This file implements the baseline sticky wage model of Jordi Galí (2010): Monetary Policy and Unemployment,
 * Handbook of Monetary Economics, Volume 3A, Chapter 10, pp. 487-546
 *
 * It demonstrates how in a linearized model a steady_state-file can be used to set the deep parameters of the
 * model to satisfy calibration targets on the non-linear model. The steady_state-file takes the calibration targets 
 * and calls a numerical solver on some of the nonlinear steady state equations to get the corresponding parameters 
 * that make the steady state satisfy the targets.
 *
 * Notes:
 * We were not able to exactly replicate the figures from the published paper. When writing this file, we encountered
 * the following issues:
 *  1. p. 516 states that Theta=0.0014. But Theta is chosen to satisfy the calibration targets given
 *          - labor share S_n=1/3
 *          - share of hiring costs to wage (W_div_PG)^-1=0.045
 *      Thus, Theta needs to satisfy
 *          Theta=delta*(W_div_PG)^-1*S_n=0.12*0.045*2/3=0.0036
 *      As the present file works with the calibration targets, the number presumably differs from Gali's
 *  2. On page 510, Upsilon is defined as the coefficient in front of the MRS and is used that 
 *      way when defining the slope of the New Keynesian Phillips Curve lambda_w. But in the Appendix on p. 542, 
 *      (1-Upsilon) instead of Upsilon denotes the coefficient on the MRS. For consistency, the formulation in 
 *      the main text is used.
 *  3. The described numbers in Gali (2010) for psi and chi do not satisfy the nonlinear steady state conditions (48-52).
 *      The present mod-file solves these equations for the correct parameters using a numerical solver
 *  4. Equation (38) should have big Phi instead of phi after (1-Upsilon)
 *  5. The definition of the unemployment rate on p. 541 
 *          urhat=fhat-nhat;
 *      implies that UR=F/N. But it should be UR=U/F, i.e. 
 *          urhat=uhat-fhat
 *      which is used here
 *  6. p. 516: in the definition of delta, it should be x/(1-x)*U/N, i.e. the brackets are missing
 * 
 * THIS MOD-FILE REQUIRES DYNARE 4.5 (I.E. THE CURRENT UNSTABLE VERSION)
 *
 * This implementation was written by Lahcen Bounader and Johannes Pfeifer. In case you spot mistakes,
 * email Johannes Pfeifer at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2016 Lahcen Bounader and Johannes Pfeifer
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


%--------------------------------------------------------------------------
                            % Endogenous Varibales
%--------------------------------------------------------------------------
@#define low_psi_calibration=0

var y_gap       ${\hat y}$              (long_name='output')
    chat        ${\hat c}$              (long_name='consumption')
    rhat        ${\hat r}$              (long_name='real interest rate')
    ihat        ${\hat i}$              (long_name='nominal interest rate')
    nhat        ${\hat n}$              (long_name='employment')
    lhat        ${\hat l}$              (long_name='labor effort')
    fhat        ${\hat f}$              (long_name='labor force')
    uhat        ${\hat u}$              (long_name='unemployment')
    uhat_0      ${\hat u^0}$            (long_name='unemployment in the begining of the period')
    urhat       ${\hat {ur}}$           (long_name='unemployment rate')
    xhat        ${\hat x}$              (long_name='job finding rate')
    ghat        ${\hat g}$              (long_name='cost of hiring')
    hhat        ${\hat h}$              (long_name='new hiring')
    mu_hat      ${\hat \mu^p}$          (long_name='markup')
    hatw_real   ${\hat \omega}$         (long_name='real wage')
    bhat        ${\hat b}$              (long_name='composite auxiliary variable')
    hatw_tar    ${\hat \omega^{tar}}$   (long_name='targeted wage (wage under Nash bargaining without rigidities)')
    a           ${a}$                   (long_name='technology shock process')
    pi_w        ${\pi^w}$               (long_name='wage inflation')
    pi_p        ${\pi^p}$               (long_name='price inflation')
    nu          ${\nu}$                 (long_name='monetary policy shock')
    urate_percentage ${U}$            (long_name='Unemployment rate in percentage points')
;

%--------------------------------------------------------------------------
                            % Exogenous Variables
%--------------------------------------------------------------------------
varexo  eps_a       ${\varepsilon_a}$   (long_name='technology shock')
        eps_nu      ${\varepsilon_\nu}$ (long_name='monetary policy shock')
       
;
%--------------------------------------------------------------------------
                            % Declaration of Parameters
%--------------------------------------------------------------------------
parameters

    alfa            ${\alpha}$              (long_name='exponent of labor in the production function')
    delta           ${\delta}$              (long_name='separation rate')
    gammma          ${\gamma}$              (long_name='coefficient of hiring cost function')
    psi             ${\psi}$                (long_name='coefficient of unemployment in the labor market effort')
    betta           ${\beta}$               (long_name='discount rate')
    rho_a           ${\rho_a}$              (long_name='autocorrelation technology shock')
    rho_nu          ${\rho_\nu}$            (long_name='autocorrelation technology shock')
    varphi          ${\varphi}$             (long_name='frish elasticity of labor effort')
    xi              ${\xi}$                 (long_name='bargaining power of workers/firms coeff')
    theta_w         ${\theta_w}$            (long_name='wage rigidities')
    theta_p         ${\theta_p}$            (long_name='price rigidities')
    phi_pi          ${\phi_{\pi}}$          (long_name='Taylor rule coeff of inflation')
    phi_y           ${\phi_y}$              (long_name='Taylor rule coeff of output gap')
    Gamma           ${\Gamma}$              (long_name='proportionality coefficient hiring cost, p. 499 bottom')
    Theta           ${\Theta}$              (long_name='share of hiring costs to GDP, p.516')
    chi             ${\chi}$                (long_name='labor disutility parameter')
    N               ${N}$                   (long_name='employment rate')
    U               ${U}$                   (long_name='unemployment rate')
    F               ${F}$                   (long_name='definition labor force, p. 516')
    L               ${L}$                   (long_name='labor in utility function, eq. (52)')
    x               ${x}$                   (long_name='steady state job finding rate')
    PG_div_W        ${\frac{PG}{W}}$        (long_name='inverse of hiring cost to wage')
    S_n             ${S^n}$                 (long_name='labor income share')
        ;

%--------------------------------------------------------------------------
                        %PARAMETRIZATION
%--------------------------------------------------------------------------

N=0.59;         %p. 515
U=0.03;         %p. 515
F=0; 
x=0.7;          %set in steady state file to be consistent with N and U
alfa=1/3;       %p. 515
betta=0.99;     %p. 515   
varphi=5;       %p. 515
theta_w=0.75;   %p. 515    
theta_p=0.75;   %p. 515
gammma=1;       %p. 515
PG_div_W=0.045; %p. 515
S_n=2/3;        %p. 515
@#if low_psi_calibration==0
    xi=0.5;      %p. 515
    psi=0;       %set in steady state file to satisfy calibration target
    chi=0;       %set in steady state file to satisfy calibration target
@#else
    xi=0.05;     %p. 515
    psi=0;       %set in steady state file to satisfy calibration target
    chi=0;       %set in steady state file to satisfy calibration target
@#endif
phi_pi=1.5;     %p. 516
phi_y=0.5/4;    %p. 516
rho_a=0.9;      %p. 517
rho_nu=0.5;     %p. 517

%--------------------------------------------------------------------------
                             %MODEL
%--------------------------------------------------------------------------
model;

//Coefficient optimal participation condition, bottom p. 541 
#Xi=(xi/(PG_div_W*(1-xi)))*(theta_w/((1-theta_w)*(1-betta*theta_w*(1-delta))));

//Slope of PC, p. 498 below equation (9)
#lambda_p=(1-theta_p)*(1-betta*theta_p)/theta_p;

//equation (24) in steady state: xi*psi*MRS=(1-xi)*x/(1-x)*G implies MRS=(1-xi)/(xi*psi)*x/(1-x)*G
//p. 510 below (35): Upsilon=xi*MRS/(W/P) with MRS from previous line: Upsilon=xi*((1-xi)/(xi*psi)*x/(1-x)*G)/(W/P)=Upsilon=(1-xi)/psi*x/(1-x)*(W/(PG))^(-1)
#Upsilon=(1-xi)/psi*x/(1-x)*PG_div_W;

//p. 501: in steady state: B=(1-(1-delta)*betta)*G
//p. 502 using B: Phi=B/(W/P+B)=1/((W/P)/B+1)=1/((W/P)/((1-(1-delta)*betta)*G)+1)
#Phi=1/(PG_div_W^(-1)/((1-(1-delta)*betta))+1);

//Slope of Wage PC, p. 511 below equation (40)
#lambda_w=(1-betta*(1-delta)*theta_w)*(1-theta_w)/(theta_w*(1-(1-Upsilon)*(1-Phi)));

[name='1. Goods Market Clearing Equations']
y_gap=(1-Theta)*chat+Theta*(ghat+hhat);

[name='2. Aggregate production function']
y_gap=a+(1-alfa)*nhat;

[name='3. Aggregate hiring and employment']
delta*hhat=nhat-(1-delta)*nhat(-1);

[name='4. Hiring Cost']
ghat=gammma*xhat;

[name='5. Job finding rate']
xhat=hhat-uhat_0;

[name='6. Effective market effort']
lhat=(N/L)*nhat+(psi*U/L)*uhat;

[name='7. Labor Force']
fhat=(N/F)*nhat+(U/F)*uhat;

[name='8. Unemployment']
uhat=uhat_0-(x/(1-x))*xhat;

[name='9. Unemployment rate'] 
urhat=uhat-fhat;

[name='10. Euler equation']
chat=chat(+1)-rhat;

[name='11. Fisherian equation']
rhat=ihat-pi_p(+1);

[name='12. Inflation equation']
pi_p=betta*pi_p(+1)-lambda_p*mu_hat;

[name='13. Optimal Hiring Condition 1']
alfa*nhat=a-((1-Phi)*hatw_real+Phi*bhat)-mu_hat;

[name='14. defintion of bhat']
bhat=(1/(1-betta*(1-delta)))*ghat-(betta*(1-delta)/(1-betta*(1-delta)))*(ghat(+1)-rhat);

[name='15. Optimal participation condition']
chat+varphi*lhat=(1/(1-x))*xhat+ghat-Xi*pi_w;

[name='16. Evolution real wage']
hatw_real=hatw_real(-1)+pi_w-pi_p;

[name='17. Wage Phillips Curve']
pi_w=betta*(1-delta)*pi_w(+1)-lambda_w*(hatw_real-hatw_tar);

[name='18. Target inflation']
hatw_tar=Upsilon*(chat+varphi*lhat)+(1-Upsilon)*(-mu_hat+a-alfa*nhat);

[name='19. Interest rate rule']
ihat=phi_pi*pi_p+phi_y*y_gap+nu;

[name='19. Monetary policy shock']
nu=rho_nu*nu(-1)+eps_nu;

[name='20. Definition of technology process']
a=rho_a*a(-1)+eps_a;

urate_percentage=U/F*urhat;
% 21. 22. Efficiency Conditions']
%log(a)-alfa*nhat=(1-Omega)*(chat+varphi*lhat)+Omega*bhat;
%chat+varphi*lhat=(1/(1-0.7))*xhat+ghat;

end;

resid(1);
check;

shocks;
    var eps_nu=0.25^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    var eps_a= 1^2;
end;

write_latex_parameter_table;
write_latex_dynamic_model;
write_latex_definitions;
collect_latex_files;

%----------------------------------------------------------------
% generate IRFs for monetary policy shock, replicates Figures 2A+B/3A+B
%----------------------------------------------------------------
stoch_simul(order = 1,irf=12) y_gap urhat urate_percentage nhat fhat pi_p hatw_real;
  
