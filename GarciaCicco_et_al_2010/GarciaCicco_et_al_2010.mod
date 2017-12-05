/*
 * This file replicates the model studied in:
 * García-Cicco, Javier and Pancrazi, Roberto and Uribe, Martín (2010): "Real Business Cycles
 * in Emerging Countries", American Economic Review, 100(5), pp. 2510-2531.
 * 
 * It provides a replication code for the main results of the original paper  
 * for the case of Argentina. 
 *
 * This mod-file shows how to use the loglinear and logdata options of Dynare (implemented since Dynare 4.5).
 * 
 * Notes:
 * - The estimation results reported in Table 3 of the paper are not easily replicable. The authors are currently working on an Erratum. The standard deviations of
 *      the estimated measurement error reported in Table 3 of the paper are actually variances. This mod-file reflects this difference to the publised version
 *      by taking the square root. Additionally, the Hessian at the mode is not well-behaved, because the mode is at a corner solution. According to communications with 
 *      GPU, they used a  positive definite approximation to the non-positive inverse Hessian to make the MCMC run. Dynare will in general run into the same problem, unless
 *      mode_compute=6 is used (which is the default in this mod-file). These problems with the Hessian may result in poor convergence of the MCMC, implying that a long chain 
 *      needs to be used (or alternatively, the use_tarb option). However, the main results of the paper are unaffected by these issues.
 * - These problems only affect the estimation. The simulations are fine. Results of stoch_simul have been checked with the moments reported in the replication
 *      code available at the AER homepage
 * - The data are taken from the AER homepage and transformed into log-differences to conform to the observables.
 * - Estimation with 2 million draws takes 6-8 hours, including mode-finding.
 * - Results for most tables can be found in the oo_ structure as documented in the Dynare manual. For example:
 *       Table 4 - Second Moments can be found in oo_.PosteriorTheoreticalMoments.dsge.correlation and covariance
 *       Table 5 - Variance Decomposition can be found in oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition
 * - The definition of the trade balance in the replication codes of the original paper is slightly incorrect. The adjustment cost 
 *       term there is "PHI/2 * (kp/k*g -G)^2" but it should be "PHI/2 * (kp/k*g -G)^2*k". This error has been corrected here. As the 
 *       wrong term is 0 up to first order, it does not affect any of the results of the paper.
 * - The GPU paper states that the interest-elastic debt premium (between equations (3) and (4)) depends on $D_{t+1}$, i.e. the debt
 *      level decided upon today (as it is predetermined). However, both the replication code and the Appendix to the paper
 *      use $D_{t}$, i.e. the debt level decided upon yesterday. This replication file follows the Appendix and the replication code.
 *    
 * This implementation was written by Johannes Pfeifer, based on the replication code by Martín Uribe. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2014-17 Johannes Pfeifer
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



@#define RBC =0
//set to 1 for RBC model and to 0 for Financial Frictions Model

var c       $c$ 
    k       $k$ 
    a       $a$ 
    h       $h$ 
    d       $d$ 
    y       $y$ 
    invest  $i$  
    tb      $tb$ 
    mu_c    ${MU_C}$ 
    tb_y    ${\frac{TB}{Y}}$ 
    g_y     ${\Delta Y}$
    g_c     ${\Delta C}$
    g_invest ${\Delta I}$
    g       ${g}$
    r       ${r}$
    mu      ${\mu}$
    nu      ${\nu}$
    @#if RBC == 0
    s       ${s}$
    @# endif
; 

predetermined_variables k d;

%Define parameters
parameters beta     ${\beta}$ 
        gamma       ${\gamma}$ 
        delta       ${\delta}$
        alpha       ${\alpha}$
        psi         ${\psi}$
        omega       ${\omega}$
        theta       ${\theta}$
        phi         ${\phi}$
        dbar        ${\bar d}$
        gbar        ${\bar g}$
        rho_a       ${\rho_a}$
        rho_g       ${\rho_g}$
        rho_nu      ${\rho_\nu}$
        rho_mu      ${\rho_\mu}$
        rho_s       ${\rho_s}$
    @#if RBC == 0
        s_share     ${sshare}$
        S           ${S}$
    @# endif
;

varexo eps_a ${\varepsilon_a}$
        eps_g ${\varepsilon_g}$ 
        eps_nu ${\varepsilon_\nu}$
        eps_mu ${\varepsilon_\mu}$
    @#if RBC == 0
        eps_s ${\varepsilon_s}$
    @# endif
;

        
@#if RBC == 1
    gbar  = 1.0050; %Gross growth rate of output
    rho_g = 0.8280; %Serial correlation of innovation in permanent technology shock
    rho_a = 0.7650; %Serial correlation of transitory technology shock
    phi   = 3.3000; %Adjustment cost parameter
@# else
    gbar  = 1.009890776104921;
    rho_g = 0.323027844166870;
    rho_a = 0.864571930755821;
    phi   = 4.810804146604144;
@# endif
            
%parameters only used for Financial frictions model, irrelevant for RBC where volatilities are 0        
rho_nu =    0.850328786147732;
rho_s  = 0.205034667802314;
rho_mu = 0.906802888826967;

%From Table 2, except for psi, which is estimated for Financial Frictions model

gamma = 2; %intertemporal elasticity of substitution
delta = 1.03^4-1;%0.03; %Depreciation rate
alpha = 0.32; %Capital elasticity of the production function
omega = 1.6; %exponent of labor in utility function
theta = 1.4*omega;
beta = 0.98^4;%0.98;%discount factor
dbar = 0.007;

@#if RBC == 1
    psi = 0.001;
@# else
    psi = 2.867166241970346; %parameter governing the debt elasticity of the interest rate. 
    s_share = 0.10; %Share of public spending in GDP
@# endif
        
        
model;
#RSTAR = 1/beta * gbar^gamma; %World interest rate
%1. Interest Rate
r = RSTAR + psi*(exp(d-dbar) - 1)+exp(mu-1)-1;

%2. Marginal utility of consumption
mu_c = nu * (c - theta/omega*h^omega)^(-gamma);

%3. Resource constraint (see the remark on the fixed typo in the preamble)
@#if RBC == 1
    y= log(tb) + c + invest + phi/2 * (k(+1)/k*g -gbar)^2*k;
@# else
    y= log(tb) + c + s + invest + phi/2 * (k(+1)/k*g -gbar)^2*k;
@# endif

%4. Trade balance
log(tb)= d - d(+1)*g/r;

%5. Definition output
y= a*k^alpha*(g*h)^(1-alpha);

%6. Definition investment
invest= k(+1)*g - (1-delta) *k;

%7. Euler equation
mu_c= beta/g^gamma*r*mu_c(+1);

%8. First order condition labor
theta*h^(omega-1)=(1-alpha)*a*g^(1-alpha)*(k/h)^alpha;

%9. First order condition investment
mu_c*(1+phi*(k(+1)/k*g-gbar))= beta/g^gamma*mu_c(+1)*(1-delta+alpha*a(+1)*(g(+1)*h(+1)/k(+1))^(1-alpha) +phi*k(+2)/k(+1)*g(+1)*(k(+2)/k(+1)*g(+1)-gbar) - phi/2*(k(+2)/k(+1)*g(+1)-gbar)^2);

%10. Definition trade-balance to output ratio
log(tb_y) = log(tb)/y; 

%11. Output growth
g_y= y/y(-1)*g(-1);

%12. Consumption growth
g_c = c/c(-1)*g(-1);

%13. Investment growth
g_invest = invest/invest(-1)*g(-1);

%14. LOM temporary TFP
log(a)=rho_a * log(a(-1))+eps_a; 

%15. LOM permanent TFP Growth
log(g/gbar)=rho_g*log(g(-1)/gbar)+eps_g; 

%16. Preference shock
log(nu) =rho_nu * log(nu(-1))+eps_nu;

%17. exogenous stochastic country premium shock
log(mu)= rho_mu * log(mu(-1))+eps_mu;

@#if RBC == 1
    
@# else
    %18. Exogenous spending shock
    log(s/S)= rho_s * log(s(-1)/S) + eps_s;
@# endif

end;

steady_state_model;
    r   = 1/beta*gbar^gamma; %World interest rate
    d   = dbar; %foreign debt
    k_over_gh =((gbar^gamma/beta-1+delta)/alpha)^(1/(alpha-1)); %k/(g*h)
    h   = ((1-alpha)*gbar*k_over_gh^alpha/theta)^(1/(omega-1)); %hours
    k   = k_over_gh*gbar*h; %capital
    invest = (gbar-1+delta)*k; %investment
    y   = k^alpha*(h*gbar)^(1-alpha); %output
    @#if RBC == 1
        s = 0;    
    @# else
        s = y*s_share;
        S = s;
    @# endif
    c   = (gbar/r-1)*d +y-s-invest; %Consumption
    tb  = y - c - s - invest; %Trade balance
    tb_y = tb /y;
    mu_c = (c - theta/omega*h^omega)^(-gamma); %marginal utility of wealth
    a   = 1; %productivity shock 
    g   = gbar; %Growth rate of nonstationary productivity shock
    g_c = g;
    g_invest = g;
    g_y = g;
    nu  = 1;
    mu  = 1;
    tb  = exp(tb);
    tb_y= exp(tb_y);
end;

shocks;
@#if RBC == 1
    var eps_a; stderr 0.0270;
    var eps_g; stderr 0.0300;
    var eps_nu; stderr  0;
    var eps_mu; stderr  0;
@# else
    var eps_a; stderr 0.033055089525252;
    var eps_g; stderr 0.010561526060797;
    var eps_nu; stderr  0.539099453618175;
    var eps_s; stderr   0.018834174505537;
    var eps_mu; stderr  0.057195449717680;
@# endif
end;

write_latex_dynamic_model; //write equations to TeX-file

stoch_simul(loglinear,order=1,irf=0) g_y g_c g_invest tb_y;

// Replicate output of replication code 
fprintf('%30s \t %5s   \t %5s   \t %5s   \t %5s   \n',' ','g_y','g_c','g_inv','TB/Y')
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Standard Deviations:',sqrt(diag(oo_.var))*100)

fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Correlation with g_y:',oo_.gamma_y{1,1}(strmatch('g_y',var_list_,'exact'),:)./sqrt(oo_.gamma_y{1,1}(strmatch('g_y',var_list_,'exact'),strmatch('g_y',var_list_,'exact')))./sqrt(diag(oo_.gamma_y{1,1}))')
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','First Order Autocorr.:',diag(oo_.autocorr{1,1}))
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Second Order Autocorr.:',diag(oo_.autocorr{1,2}))
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Third Order Autocorr.:',diag(oo_.autocorr{1,3}))
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Fourth Order Autocorr.:',diag(oo_.autocorr{1,4}))
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Correlation with TB/Y:',oo_.gamma_y{1,1}(strmatch('tb_y',var_list_,'exact'),:)./sqrt(oo_.gamma_y{1,1}(strmatch('tb_y',var_list_,'exact'),strmatch('tb_y',var_list_,'exact')))./sqrt(diag(oo_.gamma_y{1,1}))')

// Generate part of Figure 4 
verbatim;
    [data_mat,data_header]=xlsread('data_argentina.xls',1,'G2:J107');
    %sqrt(0.06*var(data_mat)); prior bounds
    figure('Name','Figure 4: Autocorrelation Function')
    tby_data=data_mat(:,strcmp('tb_y',data_header));
    plot((1:4),[corr(tby_data(2:end-3),tby_data(1:end-4)),corr(tby_data(3:end-2),tby_data(1:end-4)),corr(tby_data(4:end-1),tby_data(1:end-4)),corr(tby_data(5:end),tby_data(1:end-4))],'b-')
    hold on
    tb_pos=strmatch('tb_y',var_list_,'exact');
    plot((1:4),[oo_.autocorr{1,1}(tb_pos,tb_pos) oo_.autocorr{1,2}(tb_pos,tb_pos) oo_.autocorr{1,3}(tb_pos,tb_pos) oo_.autocorr{1,4}(tb_pos,tb_pos)],'r-.')
    xlabel('Lags')
    legend('Data','Model')
end;

% Do Table 5 with parameter of calibrated model    
fprintf('%30s \t %5s   \t %5s   \t %5s   \t %5s   \n',' ','g_y','g_c','g_inv','TB/Y')
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Non-stationary TFP:',[oo_.variance_decomposition(strmatch('g_y',var_list_,'exact'),strmatch('eps_g',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_c',var_list_,'exact'),strmatch('eps_g',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_invest',var_list_,'exact'),strmatch('eps_g',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('tb_y',var_list_,'exact'),strmatch('eps_g',M_.exo_names,'exact'))])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Stationary TFP:',[oo_.variance_decomposition(strmatch('g_y',var_list_,'exact'),strmatch('eps_a',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_c',var_list_,'exact'),strmatch('eps_a',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_invest',var_list_,'exact'),strmatch('eps_a',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('tb_y',var_list_,'exact'),strmatch('eps_a',M_.exo_names,'exact'))])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Preference:',[oo_.variance_decomposition(strmatch('g_y',var_list_,'exact'),strmatch('eps_nu',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_c',var_list_,'exact'),strmatch('eps_nu',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_invest',var_list_,'exact'),strmatch('eps_nu',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('tb_y',var_list_,'exact'),strmatch('eps_nu',M_.exo_names,'exact'))])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Country Premium:',[oo_.variance_decomposition(strmatch('g_y',var_list_,'exact'),strmatch('eps_mu',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_c',var_list_,'exact'),strmatch('eps_mu',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_invest',var_list_,'exact'),strmatch('eps_mu',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('tb_y',var_list_,'exact'),strmatch('eps_mu',M_.exo_names,'exact'))])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Ex. Spending:',[oo_.variance_decomposition(strmatch('g_y',var_list_,'exact'),strmatch('eps_s',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_c',var_list_,'exact'),strmatch('eps_s',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('g_invest',var_list_,'exact'),strmatch('eps_s',M_.exo_names,'exact')),oo_.variance_decomposition(strmatch('tb_y',var_list_,'exact'),strmatch('eps_s',M_.exo_names,'exact'))])    
    

varobs g_y g_c g_invest tb_y;
        
estimated_params;
    gbar, , , ,uniform_pdf,(1+1.03)/2,sqrt(12)^(-1)*(1.03-1),1,1.03;
    stderr eps_g, , , ,uniform_pdf,(0+0.2)/2,sqrt(12)^(-1)*(0.2-0),0,0.2;
    rho_g, , , ,uniform_pdf,(-0.99+0.99)/2,sqrt(12)^(-1)*(-0.99-0.99),-0.99,0.99;
    stderr eps_a, , , ,uniform_pdf,(0+0.2)/2,sqrt(12)^(-1)*(0.2-0),0,0.2;
    rho_a, , , ,uniform_pdf,(-0.99+0.99)/2,sqrt(12)^(-1)*(-0.99-0.99),-0.99,0.99;
    phi, , , ,uniform_pdf, (0+8)/2, sqrt(12)^(-1)*(0-8), 0, 8;

    @#if RBC == 0
        stderr eps_nu, , , ,uniform_pdf,(0+1)/2,sqrt(12)^(-1)*(1-0),0,1; //higher upper bound than the others
        rho_nu, , , ,uniform_pdf,(-0.99+0.99)/2,sqrt(12)^(-1)*(-0.99-0.99),-0.99,0.99;
        stderr eps_s, , , ,uniform_pdf,(0+0.2)/2,sqrt(12)^(-1)*(0.2-0),0,0.2;
        rho_s, , , ,uniform_pdf,(-0.99+0.99)/2,sqrt(12)^(-1)*(-0.99-0.99),-0.99,0.99;
        stderr eps_mu, , , ,uniform_pdf,(0+0.2)/2,sqrt(12)^(-1)*(0.2-0),0,0.2;
        rho_mu, , , ,uniform_pdf,(-0.99+0.99)/2,sqrt(12)^(-1)*(-0.99-0.99),-0.99,0.99;

        psi, , , ,uniform_pdf, (0+5)/2, sqrt(12)^(-1)*(0-5), 0, 5;
    @# endif

    stderr g_y, , , ,uniform_pdf, (sqrt(0.0001)+sqrt(0.13))/2,sqrt(12)^(-1)*(sqrt(0.13)-sqrt(0.0001)), sqrt(0.0001), sqrt(0.13);
    stderr g_c, , , ,uniform_pdf, (sqrt(0.0001)+sqrt(0.19))/2,sqrt(12)^(-1)*(sqrt(0.19)-sqrt(0.0001)), sqrt(0.0001), sqrt(0.19);
    stderr g_invest, , , ,uniform_pdf, (sqrt(0.0001)+sqrt(0.51))/2,sqrt(12)^(-1)*(sqrt(0.51)-sqrt(0.0001)), sqrt(0.0001), sqrt(0.51);
    stderr tb_y, , , ,uniform_pdf, (sqrt(0.0001)+sqrt(0.13))/2,sqrt(12)^(-1)*(sqrt(0.13)-sqrt(0.0001)), sqrt(0.0001), sqrt(0.13);
end;
    
estimated_params_init(use_calibration); //Use their posterior as starting values for estimation; for measurement error, only the 
    @#if RBC == 0
        stderr g_y, sqrt(0.000114861607534);
        stderr g_c,  sqrt(0.000130460798135);
        stderr g_invest,  sqrt(0.001436553959841);
        stderr tb_y,  sqrt(0.000109652068259);
    @# else
        stderr g_y,  sqrt(0.00011);
        stderr g_c,  sqrt(0.00011);
        stderr g_invest,  sqrt(0.0216);
        stderr tb_y,  sqrt(0.00011);
    @# endif
end;

estimation(datafile=data_argentina,
        xls_range=G2:J107, 
        loglinear,
        logdata, //data is already logged, loglinear option would otherwise log the data
        mode_check,
        mode_compute=6,
        moments_varendo,
        mh_nblocks=1,
        mh_replic=100000,
        consider_only_observed
                );
%Plot some parameter draws to visually check how MCMC behaved        
% trace_plot(options_,M_,estim_params_,'DeepParameter',1,'gbar');
% trace_plot(options_,M_,estim_params_,'DeepParameter',1,'rho_g');
% trace_plot(options_,M_,estim_params_,'StructuralShock',1,'eps_s');
% trace_plot(options_,M_,estim_params_,'DeepParameter',1,'phi');
% trace_plot(options_,M_,estim_params_,'MeasurementError',1,'g_invest');


fprintf('%30s \t %5s   \t %5s   \t %5s   \t %5s   \n',' ','g_y','g_c','g_inv','TB/Y')
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Standard Deviations:',[sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y),sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_c.g_c),sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_invest.g_invest),sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y)]*100)
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Correlation with g_y:',[oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y)),oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_c/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_c.g_c)),oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_invest/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_invest.g_invest)),oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.tb_y/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y))])

fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','First Order Autocorr.:',[oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_y.g_y(1),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_c.g_c(1),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_invest.g_invest(1),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.tb_y.tb_y(1)])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Second Order Autocorr.:',[oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_y.g_y(2),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_c.g_c(2),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_invest.g_invest(2),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.tb_y.tb_y(2)])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Third Order Autocorr.:',[oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_y.g_y(3),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_c.g_c(3),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_invest.g_invest(3),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.tb_y.tb_y(3)])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Fourth Order Autocorr.:',[oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_y.g_y(4),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_c.g_c(4),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.g_invest.g_invest(4),oo_.PosteriorTheoreticalMoments.dsge.correlation.Mean.tb_y.tb_y(4)])
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Correlation with TB/Y:',[oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.tb_y/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_y.g_y)),oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_c.tb_y/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_c.g_c)),oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_invest.tb_y/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.g_invest.g_invest)),oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y/(sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y)*sqrt(oo_.PosteriorTheoreticalMoments.dsge.covariance.Mean.tb_y.tb_y))])


//Do Table 5 with estimated parameters
fprintf('%30s \t %5s   \t %5s   \t %5s   \t %5s   \n',' ','g_y','g_c','g_inv','TB/Y')
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Non-stationary TFP:',[oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_y.eps_g,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_c.eps_g,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_invest.eps_g,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.tb_y.eps_g]*100)
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Stationary TFP:',[oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_y.eps_a,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_c.eps_a,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_invest.eps_a,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.tb_y.eps_a]*100)
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Preference:',[oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_y.eps_nu,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_c.eps_nu,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_invest.eps_nu,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.tb_y.eps_nu]*100)
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Country Premium:',[oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_y.eps_mu,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_c.eps_mu,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_invest.eps_mu,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.tb_y.eps_mu]*100)
@#if RBC == 0
fprintf('%30s \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n','Ex. Spending:',[oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_y.eps_s,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_c.eps_s,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.g_invest.eps_s,oo_.PosteriorTheoreticalMoments.dsge.VarianceDecomposition.Mean.tb_y.eps_s]*100)    
@# endif
