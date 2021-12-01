/*
 * This file replicates the model studied in:
 * Jordi Galí and Tommaso Monacelli (2005) "Monetary Policy and Exchange Rate
 * Volatility in a Small Open Economy", Review of Economic Studies 72, 707-734.
 *
 * It provides a replication of the main results of the original paper  
 * in Section 5 (SIMPLE MONETARY POLICY RULES FOR THE SMALL OPEN ECONOMY). It replicates:
 *  - Figure 1: Impulse responses to a domestic productivity shock under alternative policy rules
 *  - Table 1: Cyclical properties of alternative policy regimes
 *  - Table 2: Contribution to welfare losses
 * 
 * Notes:
 * - The mod-file has been tested with Dynare 4.6
 * - To generate the replication Tables and the Figure, run all four policy cases by setting the respective indicators (OPTIMAL,DITR,CITR,PEG)
 *      in turn to 1. When the code detects all four saved results files, the figure and tables will be displayed
 * - The IRFs in Figure 1 were generated with a persistence of 0.9. This can be inferred from the optimal case
 *      where the terms of trade follow the technology process one to one. 
 * - Table 2 reports the welfare losses with the leading minus sign
 * - The TeX-compilation requires PDFLaTeX to be installed and in the path
 *
 * This implementation was written by Johannes Pfeifer. If you spot mistakes, email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2015 Johannes Pfeifer
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


@#define OPTIMAL =1
@#define DITR =0
@#define CITR =0
@#define PEG =0

% define string for saving different policies
@#if DITR == 1
    case_title='domestic inflation-based Taylor rule (DITR)';
@#else
    @#if CITR ==1
         case_title='CPI inflation-based Taylor rule (CITR)';
    @#else
        @#if PEG ==1
             case_title='exchange rate peg (PEG)';
        @#else
            @#if OPTIMAL ==1
                 case_title='Optimal Policy';
            @#else
                error('One case must be set to 1')
            @#endif
        @#endif
    @#endif
@#endif


var pih ${\pi_h}$ (long_name='Domestic inflation')
    x $x$ (long_name='Output gap')
    y $y$ (long_name='Output')
    ynat ${\bar y}$ (long_name='Natural output')
    rnat ${\bar r}$ (long_name='Natural interest rate')
    r $r$ (long_name='Nominal interest rate')
    s $s$ (long_name='Terms of trade')
    pi ${\pi}$ (long_name='CPI inflation')
    p $p$ (long_name='CPI level')
    ph ${p_h}$ (long_name='Domestic price level')
    e $e$ (long_name='Exchange rate')
    ystar ${y^*}$ (long_name='World output')
    pistar ${\pi^{*}}$ (long_name='World inflation')
    n ${n}$ (long_name='Employment')
    nx ${nx}$ (long_name='Net Exports')
    real_wage ${w-p}$ (long_name='Real Wage')
    a $a$ (long_name='Risk aversion')
    c $c$ (long_name='Domestic consumption')
    deprec_rate $\Delta e_t$ (long_name='Nominal depr. rate')
    ;

varexo eps_star ${\varepsilon^{*}}$ (long_name='World output shock')
    eps_a ${\varepsilon^{a}}$ (long_name='World output shock');

parameters sigma $\sigma$ (long_name='risk aversion')
    eta $\eta$ (long_name='Substitution home foreign')
    gamma $\gamma$ (long_name='Substitution between foreign')
    phi $\varphi$ (long_name='Inverse Frisch elasticity')
    epsilon $\varepsilon$ (long_name='Elasticit of substitution')
    theta $\theta$ (long_name='Calvo parameter')
    beta $\beta$ (long_name='discount factor')
    alpha $\alpha$ (long_name='openness')
    phi_pi $\phi_\pi$ (long_name='Feedback Taylor rule inflation')
    rhoa $\rho_a$ (long_name='autocorrelation TFP')
    rhoy $\rho_y$ (long_name='autocorrelation foreign output')
;

% set deep parameters
sigma = 1;
eta = 1 ;
gamma = 1;
phi = 3;
epsilon = 6;
theta = 0.75;
beta  = 0.99;
alpha = 0.4;
phi_pi = 1.5;
rhoa = 0.9; %use value used for Figure 1, reset later                                                                
rhoy = 0.86;  

model(linear);
//define parameter dependencies
//steady state real interest rate, defined below equation (11)
#rho  = beta^(-1)-1;
//defined below equation (27)
#omega = sigma*gamma+(1-alpha)*(sigma*eta-1);
//defined below equation (29)
#sigma_a =sigma/((1-alpha)+alpha*omega);

#Theta=(sigma*gamma-1)+(1-alpha)*(sigma*eta-1);
//defined below equation (32)
#lambda = (1-(beta*theta))*(1-theta)/theta;
//defined below equation (36)
#kappa_a =lambda*(sigma_a+phi);
//defined below equation (35)
#Gamma = (1+phi)/(sigma_a+phi);
#Psi = -Theta*sigma_a/(sigma_a+phi);

[name='Equation (37), IS Curve']
x    = x(+1) - sigma_a^(-1)*(r - pih(+1) - rnat) ;                              
[name='Equation (36), Philips Curve']
pih  = beta * pih(+1)+ kappa_a*x;                                                
[name='Equation below (37)']
rnat = -sigma_a*Gamma*(1-rhoa)*a + alpha*sigma_a*(Theta+Psi)*(ystar(+1)-ystar);
[name='Equation (35), definition natural level of output']
ynat = Gamma*a + alpha*Psi*ystar;                                                 
[name='Equation above (35), definition output gap']
x    = y - ynat;                                                               
[name='Equation (29)']
y = ystar + sigma_a^(-1)*s;
[name='Equation (14)']
pi   = pih + alpha*(s-s(-1));
[name='Equation 15 (first difference)']
s    = s(-1) + e - e(-1) + pistar - pih;
[name='Constant world inflation, see p.724 (Given constant world prices)']
pistar = 0;
[name='Equation (22), employment']
y = a + n;
[name='Equation (31), net exports']
nx = alpha*(omega/sigma-1)*s;
[name='Equation (27), defines consumption']
y = c+alpha*omega/sigma*s;
[name='Above equation (11), defines real wage']
real_wage = sigma*c+phi*n;

[name='stochastic process for technology p. 723']
a    = rhoa*a(-1) + eps_a;
[name='stochastic process for foreign output p. 723']
ystar= rhoy*ystar(-1) + eps_star;

// Equations on page 
@#if DITR == 1
[name='domestic inflation-based Taylor rule (DITR)']
r = phi_pi*pih;
@#else
    @#if CITR ==1
    [name='CPI inflation-based Taylor rule (CITR)']
        r = phi_pi*pi;
    @#else
        @#if PEG ==1
            [name='exchange rate peg (PEG)']
                e=0;
        @#else
            @#if OPTIMAL ==1
                [name='optimal policy']
                pih=0;
            @#else
            
            @#endif
        @#endif
    @#endif
@#endif
  
[name='definition consumer price level']
pi   = p - p(-1);
[name='definition domestic price level']
pih  = ph - ph(-1);
[name='definition nominal depreciation rate of exchange rate']
deprec_rate=e-e(-1);
end;

shocks;
var eps_a = 1; //unit shock
end;

//generate LaTeX-files with equations and parameterization
write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

stoch_simul(TeX,order=1,irf=20,irf_plot_threshold=0) pih x pi s e r ph p a;

collect_latex_files;

//Uncomment the following lines to generate a PDF file using PDFLaTeX (if installed)
% if system(['pdflatex -halt-on-error ' M_.fname '_TeX_binder.TeX'])
%     error('TeX-File did not compile.')
% end

set_param_value('rhoa',0.66); %reset to value stated on page 723 and used for                                                                

shocks;
var eps_a = 0.0071^2; //the standard deviation of productivity
var eps_star = 0.0078^2; //the standard deviation of world output
var eps_a, eps_star = 0.3*0.0071*0.0078; //covariance, constructed from correlation and standard deviations                                            
end;

stoch_simul(order=1,irf=0);


var_string={'y','pih','pi','r','s','deprec_rate'};

fprintf('\nTABLE 1: Cyclical properties of alternative policy regimes\n')
fprintf('Case: %s\n', case_title)
for var_iter=1:length(var_string)
    var_pos=strmatch(var_string{var_iter},M_.endo_names,'exact');
    cyc_moments(var_iter,1)=sqrt(oo_.var(var_pos,var_pos))*100;
    fprintf('%20s \t %3.2f \n',M_.endo_names_long{strmatch(var_string{var_iter},M_.endo_names,'exact'),:},cyc_moments(var_iter,1))
end


%%%%%%%%%%%%%%%% Do entries of Table 2 %%%%%%%%%%%%%%%%
%find output gap and inflation in covariance matrix
x_pos=strmatch('x',M_.endo_names,'exact');
pih_pos=strmatch('pih',M_.endo_names,'exact');


%%%%%%%%%%%%Panel 1
V_infl(1,1)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos))*100;
V_output_gap(1,1)=-(1-alpha)/2*((1+phi)*oo_.var(x_pos,x_pos))*100;
V(1,1)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos)+(1+phi)*oo_.var(x_pos,x_pos))*100;

fprintf('\nTABLE 2: Contribution to welfare losses\n')
fprintf('Case: %s\n', case_title)
fprintf('mu=%2.1f, phi=%2.1f\n',M_.params(strmatch('epsilon',M_.param_names,'exact'))/(M_.params(strmatch('epsilon',M_.param_names,'exact'))-1),M_.params(strmatch('phi',M_.param_names,'exact')))
fprintf('%-20s \t %5.4f \n','Var(Domestic infl.)',V_infl(1,1))
fprintf('%-20s \t %5.4f \n','Var(Output gap)',V_output_gap(1,1))
fprintf('%-20s \t %5.4f \n\n','Total',V(1,1))

%%%%%%%%%%%%Panel 2
mu=1.1;
set_param_value('epsilon',mu/(mu-1)); 
set_param_value('phi',3); 

stoch_simul(order=1,irf=0,noprint);

V_infl(1,2)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos))*100;
V_output_gap(1,2)=-(1-alpha)/2*((1+phi)*oo_.var(x_pos,x_pos))*100;
V(1,2)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos)+(1+phi)*oo_.var(x_pos,x_pos))*100;

fprintf('mu=%2.1f, phi=%2.1f\n',M_.params(strmatch('epsilon',M_.param_names,'exact'))/(M_.params(strmatch('epsilon',M_.param_names,'exact'))-1),M_.params(strmatch('phi',M_.param_names,'exact')))
fprintf('%-20s \t %5.4f \n','Var(Domestic infl.)',V_infl(1,2))
fprintf('%-20s \t %5.4f \n','Var(Output gap)',V_output_gap(1,2))
fprintf('%-20s \t %5.4f \n\n','Total',V(1,2))

%%%%%%%%%%%%Panel 3
mu=1.2;
set_param_value('epsilon',mu/(mu-1)); %value page 723    
set_param_value('phi',10); %value page 723    

stoch_simul(order=1,irf=0,noprint);
V_infl(1,3)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos))*100;
V_output_gap(1,3)=-(1-alpha)/2*((1+phi)*oo_.var(x_pos,x_pos))*100;
V(1,3)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos)+(1+phi)*oo_.var(x_pos,x_pos))*100;

fprintf('mu=%2.1f, phi=%2.1f\n',M_.params(strmatch('epsilon',M_.param_names,'exact'))/(M_.params(strmatch('epsilon',M_.param_names,'exact'))-1),M_.params(strmatch('phi',M_.param_names,'exact')))
fprintf('%-20s \t %5.4f \n','Var(Domestic infl.)',V_infl(1,3))
fprintf('%-20s \t %5.4f \n','Var(Output gap)',V_output_gap(1,3))
fprintf('%-20s \t %5.4f \n\n','Total',V(1,3))

%%%%%%%%%%%%Panel 4
mu=1.1;
set_param_value('epsilon',mu/(mu-1)); %value page 723    
set_param_value('phi',10); %value page 723    

stoch_simul(order=1,irf=0,noprint);

V_infl(1,4)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos))*100;
V_output_gap(1,4)=-(1-alpha)/2*((1+phi)*oo_.var(x_pos,x_pos))*100;
V(1,4)=-(1-alpha)/2*(epsilon/((1-(beta*theta))*(1-theta)/theta)*oo_.var(pih_pos,pih_pos)+(1+phi)*oo_.var(x_pos,x_pos))*100;

fprintf('mu=%2.1f, phi=%2.1f\n',M_.params(strmatch('epsilon',M_.param_names,'exact'))/(M_.params(strmatch('epsilon',M_.param_names,'exact'))-1),M_.params(strmatch('phi',M_.param_names,'exact')))
fprintf('%-20s \t %5.4f \n','Var(Domestic infl.)',V_infl(1,4))
fprintf('%-20s \t %5.4f \n','Var(Output gap)',V_output_gap(1,4))
fprintf('%-20s \t %5.4f \n\n','Total',V(1,4))

%Save results
@#if DITR == 1
    save('Gali_Monacelli_2005_DITR','V_infl','V_output_gap','V','cyc_moments','M_','oo_','options_','case_title');
@#else
    @#if CITR ==1
         save('Gali_Monacelli_2005_CITR','V_infl','V_output_gap','V','cyc_moments','M_','oo_','options_','case_title');
    @#else
        @#if PEG ==1
             save('Gali_Monacelli_2005_PEG','V_infl','V_output_gap','V','cyc_moments','M_','oo_','options_','case_title');
        @#else
            @#if OPTIMAL ==1
                 save('Gali_Monacelli_2005_OPTIMAL','V_infl','V_output_gap','V','cyc_moments','M_','oo_','options_','case_title');
            @#else
                error('Undefined case')
            @#endif
        @#endif
    @#endif
@#endif

%%%%%%%%%%%%%%%%%%% generate Figure and Tables when all results files exist
if exist('Gali_Monacelli_2005_DITR.mat') && exist('Gali_Monacelli_2005_CITR.mat') && exist('Gali_Monacelli_2005_PEG.mat') && exist('Gali_Monacelli_2005_OPTIMAL.mat')
    DITR=load('Gali_Monacelli_2005_DITR.mat');
    CITR=load('Gali_Monacelli_2005_CITR.mat');
    PEG=load('Gali_Monacelli_2005_PEG.mat');
    OPTIMAL=load('Gali_Monacelli_2005_OPTIMAL.mat');
    var_string={'pih','x','pi','s','e','r','ph','p'};
    figure
    for fig_iter=1:length(var_string)
        subplot(4,2,fig_iter)
        plot(1:20,OPTIMAL.oo_.irfs.([var_string{fig_iter},'_eps_a']),'b-',1:20,DITR.oo_.irfs.([var_string{fig_iter},'_eps_a']),'g--',1:20,CITR.oo_.irfs.([var_string{fig_iter},'_eps_a']),'r-x',1:20,PEG.oo_.irfs.([var_string{fig_iter},'_eps_a']),'c-s')
        grid on
        title(OPTIMAL.M_.endo_names_long(strmatch(var_string{fig_iter},OPTIMAL.M_.endo_names,'exact'),:))
    end

    var_string={'y','pih','pi','r','s','deprec_rate'};

    fprintf('\nTABLE 1: Cyclical properties of alternative policy regimes\n')
    fprintf('%20s \t %3s \t %3s \t %3s \t %3s \n','','Opt','DIT','CIT','PEG')
    fprintf('%20s \t %3s \t %3s \t %3s \t %3s \n','','sd%','sd%','sd%','sd%')
    for var_iter=1:length(var_string)
        fprintf('%20s \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n',M_.endo_names_long(strmatch(var_string{var_iter},M_.endo_names,'exact'),:),OPTIMAL.cyc_moments(var_iter),DITR.cyc_moments(var_iter),CITR.cyc_moments(var_iter),PEG.cyc_moments(var_iter))
    end
    
    case_string={'Benchmark mu = 1.2, phi = 3', 'Low steady state mark-up mu = 1.1, phi = 3','Low elasticity of labour supply mu = 1.2, phi = 10','Low mark-up and elasticity of labour supply mu = 1.1, phi= 10'};
    fprintf('\nTABLE 2: Contribution to welfare losses\n')
    fprintf('%20s \t %6s \t %6s \t %6s \n','','DIT','CIT','PEG')
    for case_iter=1:4
        fprintf('%-s \n',case_string{case_iter});
        fprintf('%-20s \t %5.4f \t %5.4f \t %5.4f \n','Var(Domestic infl.)',DITR.V_infl(1,case_iter),CITR.V_infl(1,case_iter),PEG.V_infl(1,case_iter))
        fprintf('%-20s \t %5.4f \t %5.4f \t %5.4f \n','Var(Output gap)',DITR.V_output_gap(1,case_iter),CITR.V_output_gap(1,case_iter),PEG.V_output_gap(1,case_iter))
        fprintf('%-20s \t %5.4f \t %5.4f \t %5.4f \n\n','Total',DITR.V(1,case_iter),CITR.V(1,case_iter),PEG.V(1,case_iter))
    end   
end
