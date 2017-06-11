/*
 * This file implements the baseline sticky wage model of Jordi Galí (2010): Monetary Policy and Unemployment,
 * Handbook of Monetary Economics, Volume 3A, Chapter 10, pp. 487-546
 *
 * Special thanks go to Jordi Gali for providing his original codes, which allowed to clarify important calibration questions.
 *
 * Notes:
 *  1. While the figures in the paper can be replicated, they are strictly speaking not consistent with the stated calibration
 *      targets:
 *      a) p. 516 states that Theta=0.0014. But Theta is in principle chosen to satisfy the calibration targets given
 *          - labor share S_n=2/3
 *          - share of hiring costs to wage (W_div_PG)^-1=0.045
 *         Taking the formula at face value, Theta would need to satisfy
 *          Theta=delta*N*G/Y=delta*W/P/(W/P)*N/Y*G=delta*(W/(PG))^-1*S_n=0.12*0.045*2/3=0.0036
 *         The present mod-file uses Theta=0.0014.
 *      b) The calibration treats the labor share S_n as a free parameter in Theta, while it is endogenous to the 
 *          model calibration as S_n=W/P*N/Y. But the above expansion of Theta relies on the actual labor share. 
 *          Thus, when continuing the model calibration under the false pretense that in the model W/P*N/Y is 
 *          actually equal to S_n=2/3, the resulting ratio of hiring costs to GDP in the model 
 *          is 0.018 instead of 0.045.
 *  2. On page 510, Upsilon is defined as the coefficient in front of the MRS and is used that 
 *      way when defining the slope of the New Keynesian Phillips Curve lambda_w. But in the Appendix on p. 542, 
 *      (1-Upsilon) instead of Upsilon denotes the coefficient on the MRS. For consistency, the formulation in 
 *      the main text is used.
 *  3. Equation (38) should have big Phi instead of phi after (1-Upsilon)
 *  4. The definition of the unemployment rate on p. 541 
 *          urhat=fhat-nhat;
 *      involves an approximation error. The unemployment rate is defined as UR_t=U_t/F_t=(F_t-N_t)/U_t
 *      A first order Taylor approximation and expanding the RHS variables by their steady state yields 
 *          ur_t=U/F*u_hat_t-U/F*f_hat_t
 *      where ur_t is the absolute deviation of the unemployment rate from its steady state and therefore 
 *      in percentage points. Using the definition of the loglinearized labor force we have 
 *          U/F*u_hat_t=f_hat_t-N/F*n_hat_t
 *      this can be written as
 *          ur_t=f_hat_t-N/F*n_hat_t-U/F*f_hat_t=N/F*f_hat_t-N/F*n_hat_t
 *      As N/F is very close to 1, this is approximately
 *          urhat=uhat-fhat
 *      which is the equation used in Gali. Here, we use full
 *          ur_t=U/F*u_hat_t-U/F*f_hat_t
 *  5. p. 516: in the definition of delta, it should be x/(1-x)*U/N, i.e. the brackets are missing
 * 
 * THIS MOD-FILE REQUIRES DYNARE 4.5 OR HIGHER
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
    urhat_Gali  ${\hat {ur_{gali}}}$    (long_name='unemployment rate')
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
    varphi          ${\varphi}$             (long_name='Frisch elasticity of labor effort')
    xi              ${\xi}$                 (long_name='bargaining power of workers/firms coeff')
    theta_w         ${\theta_w}$            (long_name='wage rigidities')
    theta_p         ${\theta_p}$            (long_name='price rigidities')
    phi_pi          ${\phi_{\pi}}$          (long_name='Taylor rule coeff of inflation')
    phi_y           ${\phi_y}$              (long_name='Taylor rule coeff of output gap')
    Gamma           ${\Gamma}$              (long_name='proportionality coefficient hiring cost, p. 499 bottom')
    Theta           ${\Theta}$              (long_name='share of hiring costs to GDP, p.516')
    Upsilon         ${\Upsilon}$            (long_name='Composite parameter, eq. (35), p.510')
    Phi             ${\Phi}$                (long_name='share of hiring costs to hiring costs plus wage, eq. (13)')
    Xi              ${\Xi}$                 (long_name='Coefficient optimal participation condition, bottom p. 541')
    chi             ${\chi}$                (long_name='labor disutility parameter')
    N               ${N}$                   (long_name='employment rate')
    U               ${U}$                   (long_name='unemployment rate')
    F               ${F}$                   (long_name='definition labor force, p. 516')
    L               ${L}$                   (long_name='labor in utility function, eq. (52)')
    x               ${x}$                   (long_name='steady state job finding rate')
        ;

%--------------------------------------------------------------------------
                        %PARAMETRIZATION
%--------------------------------------------------------------------------

N=0.59;         %p. 515
U=0.03;         %p. 515
F=0;            %set in steady state file to be consistent with N and U
x=0.7;          
alfa=1/3;       %p. 515
betta=0.99;     %p. 515   
varphi=5;       %p. 515
theta_w=0.75;   %p. 515    
theta_p=0.75;   %p. 515
gammma=1;       %p. 515
Theta=0.0014;   %p. 515; in principle set to satisfy calibration target, see Header 
@#if low_psi_calibration==0
    xi=0.5;      %p. 515
    psi=0;       %set in steady state file to satisfy calibration target
    chi=0;       %set in steady state file to satisfy calibration target
@#else
    xi=0.05;     %p. 515
    psi=0;       %set in steady state file to satisfy calibration target
    chi=0;       %set in steady state file to satisfy calibration target
@#endif
phi_pi=1.5;     %p. 521
phi_y=0.5/4;    %p. 516
rho_a=0.9;      %p. 517
rho_nu=0.5;     %p. 517
Upsilon=0;      %set in steady state file to satisfy calibration target
Phi=0;          %set in steady state file to satisfy calibration target
Xi=0;           %set in steady state file to satisfy calibration target

%--------------------------------------------------------------------------
                             %MODEL
%--------------------------------------------------------------------------
model;

//Slope of PC, p. 498 below equation (9)
#lambda_p=(1-theta_p)*(1-betta*theta_p)/theta_p;

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
urhat=U/F*uhat-U/F*fhat;
urhat_Gali=fhat-nhat;

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
stoch_simul(order = 1,irf=12) y_gap urhat nhat fhat pi_p hatw_real;
  
