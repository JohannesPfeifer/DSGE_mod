/*
 * This file implements the RBC model of Jermann/Quadrini (2012): "Macroeconomic effects of
 * financial shocks", American Economic Review, 102(1): 238–271.
 *
 * It allows replicating the original results by setting replication = 1. When setting replication=0 
 * it generates the results of Pfeifer (2016): "Macroeconomic effects of financial shocks", which documents
 * a mistake in the TFP-construction of JQ that requires recalibrating the model.
 *
 * Notes: 
 *  - before running this mod-file, the construct_data.m must be executed first with the respective dummy 
 *      having been set
 *  - Requires Dynare 4.5
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

@#define replication=0
        
var c ${c}$ (long_name='Consumption')
     n ${n}$ (long_name='Labor')
     w ${w}$ (long_name='Wage')
     k ${k}$ (long_name='Capital')
     R ${R}$ (long_name='Effective Gross Interest Rate')
     r ${r}$ (long_name='Gross Interest Rate')
     d ${d}$ (long_name='Dividend')
     b ${b}$ (long_name='Bond')
     mu ${\mu}$ (long_name='Lagrange Multiplier')
     v ${v}$ (long_name='Value of the Firm')
     z ${z}$ (long_name='Technology')
     xi ${\xi}$ (long_name='Financial Conditions')
     y ${y}$ (long_name='Output')
     invest ${i}$ (long_name='Investment')
     yhat ${\hat y}$ (long_name='output deviation from trend')
     chat ${\hat c}$ (long_name='consumption deviation from trend') 
     ihat ${\hat i}$ (long_name='investment deviation from trend')
     nhat ${\hat n}$ (long_name='hours deviation from trend')
     byhat ${\hat y}$ (long_name='debt-repurchase to GDP ratio deviation from trend') 
     dyhat ${\hat y}$ (long_name='dividend to GDP ratio deviation from trend')
     muhat ${\hat \mu}$ (long_name='Lagrange multiplier from trend')
     vyhat ${\hat y}$ (long_name='Equity value deviation from trend')
         ;

varexo eps_z ${\varepsilon_z}$ (long_name='Technology shock')
     eps_xi ${\varepsilon_{\xi}}$ (long_name='Financial Shock')
         ;

parameters theta ${\theta}$ (long_name='capital share')
        betta ${\beta}$ (long_name='discount factor')
        alppha ${\alpha}$ (long_name='disutility from work')
        delta ${\delta}$ (long_name='depreciation')
        tau ${\tau}$ (long_name='tax wedge')
        kappa ${\kappa}$ (long_name='equity cost')
        siggma ${\sigma}$ (long_name='risk aversion')
        sigma_z ${\sigma_z}$ (long_name='std_z')
        sigma_xi ${\sigma_xi}$ (long_name='std_xi')
        BY_ratio ${(\bar b/(1+\bar r)/\bar Y}$ (long_name='Debt output ratio')
        A11 ${A_{11}}$ (long_name='A_11')
        A12 ${A_{12}}$ (long_name='A_12')
        A21 ${A_{21}}$ (long_name='A_21')
        A22 ${A_{22}}$ (long_name='A_22')
        ;
     
siggma=1;
theta = 0.36;
betta = 0.9825;
delta = 0.025;
tau = 0.35;        
BY_ratio=3.36;


@#if replication==1
    kappa = 0.146;       
    kappa_store=kappa;
    sigma_xi = 0.0098;  
    sigma_z = 0.0045; 
    covariance_z_xi =0;
    A11 = 0.9457;
    A12 = -0.0091;
    A21 = 0.0321;
    A22 = 0.9703;
@#else    
    kappa = 0.08;       
    kappa_store=kappa;
    estimated_process=load('innovations_corrected.mat')
    cov_matrix=estimated_process.cov_matrix;
    sigma_xi = sqrt(cov_matrix(2,2));  
    sigma_z = sqrt(cov_matrix(1,1)); 
    covariance_z_xi =0; %set to 0 in JQ; otherwise use covariance_z_xi =cov_matrix(2,1);
    bhat_OLS=estimated_process.bhat_OLS;
    A11 = bhat_OLS(1,1);
    A12 = bhat_OLS(1,2);
    A21 = bhat_OLS(2,1);
    A22 = bhat_OLS(2,2);    
@#endif
options_.TeX=1;

model;
[name='FOC labor, equation 1 on p. 4 Appendix']
w/c^siggma-alppha/(1-n)=0;
[name='Euler equatoin, equation 2 on p. 4 Appendix']
c^(-siggma)=betta*((R-tau)/(1-tau))*c(+1)^(-siggma);
[name='budget constraint household, equation 3 on p. 4 Appendix']
w*n+b(-1)-b/R+d-c=0;
[name='FOC labor input, equation 4 on p. 4 Appendix']
(1-theta)*z*k(-1)^theta*n^(-theta)=w*(1/(1-mu*(1+2*kappa*(d-steady_state(d)))));
[name='FOC capital, equation 5 on p. 4 Appendix']
betta*(c/(c(+1)))^siggma*((1+2*kappa*(d-steady_state(d)))/(1+2*kappa*(d(+1)-steady_state(d))))*(1-delta+(1-mu(+1)*(1+2*kappa*(d(+1)-steady_state(d))))*theta*z(+1)*k^(theta-1)*n(+1)^(1-theta))+xi*mu*(1+2*kappa*(d-steady_state(d)))=1;
[name='FOC bonds, equation 6 on p. 4 Appendix']
R*betta*(c/(c(+1)))^siggma*((1+2*kappa*(d-steady_state(d)))/(1+2*kappa*(d(+1)-steady_state(d))))+xi*mu*(1+2*kappa*(d-steady_state(d)))*(R*(1-tau)/(R-tau))=1;
[name='budget constraint firm, equation 7 on p. 4 Appendix']
(1-delta)*k(-1)+z*k(-1)^theta*n^(1-theta)-w*n-b(-1)+b/R-k-(d+kappa*(d-steady_state(d))^2)=0;
[name='Enforcement constraint, equation 8 on p. 4 Appendix']
xi*(k-b*((1-tau)/(R-tau)))=z*k(-1)^theta*n^(1-theta);

[name='Equation 1 of VAR in equation (11)']
log(z/steady_state(z))=A11*log(z(-1)/steady_state(z))+A12*log(xi(-1)/steady_state(xi))+eps_z;
[name='Equation 2 of VAR in equation (11)']
log(xi/steady_state(xi))=A21*log(z(-1)/steady_state(z))+A22*log(xi(-1)/steady_state(xi))+eps_xi;
        
[name='Production function']
y=z*k(-1)^theta*n^(1-theta);
[name='LOM capital']
invest=k-(1-delta)*k(-1);
[name='Definition firm value, equation 11) in appendix p.']
v=d+betta*c/c(+1)*v(+1);
[name='Definition before tax interest rate']
r=(R-tau)/(1-tau)-1;
[name='Definition output deviations from trend']
yhat=log(y)-log(steady_state(y));
[name='Definition consumption deviation from trend']
chat=log(c)-log(steady_state(c));
[name='Definition investment deviation from trend']
ihat=log(invest)-log(steady_state(invest));
[name='Definition hours deviation from trend']
nhat=log(n)-log(steady_state(n));
[name='Definition lagrange multiplier deviation from trend']
muhat=log(mu)-log(steady_state(mu));
[name='Definition debt repurchase share in GDP']
byhat=(b(-1)/(1+r(-1))-b/(1+r))/y;            
[name='Definition equity payout to GDP ratio']
dyhat=d/y;
[name='Definition equity share']
vyhat=log(v/(k(-1)-b(-1)))-log(steady_state(v/(k(-1)-b(-1))));
end;

shocks;
    var eps_xi=sigma_xi^2;
    var eps_z= sigma_z^2;       
    var eps_z,eps_xi = covariance_z_xi;
end;

steady;

write_latex_original_model;
write_latex_dynamic_model;
write_latex_parameter_table;

%% make VAR diagonal for generation of IRFS
        
set_param_value('A12',0);
set_param_value('A21',0);

stoch_simul(order=1,bandpass_filter=[6,32],irf=50,nograph,nocorr,nofunctions);

%% Display some statistics
fprintf('(b/(1+r)/Y: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))/(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact'))))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))
fprintf('Debt-Capital Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))/(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact'))))/oo_.dr.ys(strmatch('k',M_.endo_names,'exact')))
fprintf('Total Debt-Output Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))+oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))/(oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))*(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact')))))
fprintf('Total Debt-Capital Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))+oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))/(oo_.dr.ys(strmatch('k',M_.endo_names,'exact'))*(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact')))))
fprintf('Capital-Output Ratio: %4.3f\n',(oo_.dr.ys(strmatch('k',M_.endo_names,'exact'))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))))


%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Figure 6: Impulse Responses')
        subplot(2,4,1)
        plot(-oo_.irfs.yhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.yhat_eps_xi*100,'r--')
        axis([-inf inf -0.8 0]);
        ll=legend('TFP shock','Financial Shock');
        title('Output')
        
        subplot(2,4,2)
        plot(-oo_.irfs.nhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.nhat_eps_xi*100,'r--')
        axis([-inf inf -1.2 0.8]);
        title('Hours')
        
        subplot(2,4,3)
        plot(-oo_.irfs.chat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.chat_eps_xi*100,'r--')
        axis([-inf inf -0.6 0.06]);
        title('Consumption')
        
        subplot(2,4,4)
        plot(-oo_.irfs.ihat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.ihat_eps_xi*100,'r--')
        axis([-inf inf -3.5 0.5]);
        title('Investment')
        
        subplot(2,4,5)
        plot(-oo_.irfs.byhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.byhat_eps_xi*100,'r--')
        axis([-inf inf -1.9 2.4]);
        title('Debt rep./Y')
        
        subplot(2,4,6)
        plot(-oo_.irfs.dyhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.dyhat_eps_xi*100,'r--')
        axis([-inf inf -1.5 1.5]);
        title('Equity Payout/Y')
        
        subplot(2,4,7)
        plot(-oo_.irfs.vyhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.vyhat_eps_xi*100,'r--')
        axis([-inf inf -0.3 0.4]);
        title('Equity Value')
        
        subplot(2,4,8)
        plot(-oo_.irfs.muhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.muhat_eps_xi*100,'r--')
        axis([-inf inf -20 25])
        title('Multiplier, \mu')

%%%%%%%%%%%%%%%%%%%%%%%% Do simulations for counterfactuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#if replication==1
    innovations=load('innovations_replication.mat','resid');           
    emp_data=load('emp_data.mat');                                 % storing empirical data 
@#else
    innovations=load('innovations_corrected.mat');
    emp_data=load('emp_data.mat');                                 % storing empirical data 
@#endif

%% empirical data 
emp_GDP=emp_data.log_Real_GDP_detrended;    
emp_value_added=emp_data.log_Real_Business_value_added_detrended(:,1);
emp_equity=emp_data.equity_payout_detrended;
emp_debt=emp_data.debt_repurchases_detrended;
emp_hours=emp_data.log_Total_Private_hours_detrended;
emp_consumption=emp_data.log_Real_Personal_Consumption;
emp_investment=emp_data.log_Real_Investment_detrended;
timeline=emp_data.data_timeline;  
n_points=length(emp_GDP);      

%% Do "no frictions"-case first
        
set_param_value('kappa',0);
set_param_value('tau',0.00001);        % trick to leave mu determined in steady state (if set to zero, log of steady state of muhat (log(0)) is NaN)
@#if replication==1
    set_param_value('A12',-0.0091);
    set_param_value('A21',0.0321);
@#else
    set_param_value('A12',estimated_process.bhat_OLS(1,2));
    set_param_value('A21',estimated_process.bhat_OLS(2,1));
@#endif
                                        
stoch_simul(order=1,irf=105,nograph,nomoments,nocorr,nofunctions);

%% initialize IRF generation
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag);
shock_matrix = zeros(n_points,M_.exo_nbr);      %create shock matrix with number of time periods in colums

%% set shocks for 'productivity shocks only' 
        
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = zeros(1,n_points); 
y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_prod_nofric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);     % deviation from steady state  
        
%% set shocks for 'financial shocks only' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = zeros(1,n_points);
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_fin_nofric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);         
        
%% set shocks for 'both shocks' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_both_nofric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);


%% Now do model with frictions (baseline model)

@#if replication==1
    set_param_value('kappa',kappa_store);
    set_param_value('tau',0.35);
@#else    
    set_param_value('kappa',kappa_store);
    set_param_value('tau',0.35);
@#endif

stoch_simul(order=1,bandpass_filter=[6,32],irf=105,nograph,nomoments,nocorr,nofunctions);

%% initialize IRF generation
initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag);
shock_matrix = zeros(n_points,M_.exo_nbr);                
        
%% set shocks for 'productivity shocks only' 

shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = zeros(1,n_points); 
y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_prod_fric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);               
        
%% set shocks for 'technology shocks only' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = zeros(1,n_points);
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_fin_fric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);                    
        
%% set shocks for 'both shocks' 
shock_matrix(:,strmatch('eps_z',M_.exo_names,'exact')) = innovations.resid(:,1)';
shock_matrix(:,strmatch('eps_xi',M_.exo_names,'exact')) = innovations.resid(:,2)'; 
y2 = simult_(initial_condition_states,oo_.dr,shock_matrix,1);
y_IRF_both_fric = y2(:,M_.maximum_lag+1:end)-repmat(oo_.dr.ys,1,n_points);     
        
        
%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Figure 2: Multiplier')

plot(timeline,100*y_IRF_both_fric(strmatch('muhat',M_.endo_names,'exact'),:));
axis([-inf inf -90 120]);
title('Lagrange Multiplier');

%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff=figure('Name','Figure 3: Response to productivity shock only');

subplot(2,2,1)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('yhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*y_IRF_prod_nofric(strmatch('yhat',M_.endo_names,'exact'),:),'r--');
hh3=plot(timeline,100*emp_GDP,'g');
hold off
axis([-inf inf -14 8]);
plot_NBER_recessions([hh1;hh2;hh3]);     %plot shaded recession dates behind plot with handle h
ll=legend([hh1 hh2 hh3],'Baseline','No fin. fric.','Data');
set(ll,'Location','NorthWest');
title('GDP');

subplot(2,2,2)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('nhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*y_IRF_prod_nofric(strmatch('nhat',M_.endo_names,'exact'),:),'r--');
hh3=plot(timeline,100*emp_hours,'g');
hold off
axis([-inf inf -14 8]);        
plot_NBER_recessions([hh1;hh2;hh3]);     %plot shaded recession dates behind plot with handle h
title('Hours');

subplot(2,2,3)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('byhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_debt,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Debt repurchase');

subplot(2,2,4)
hh1=plot(timeline,100*y_IRF_prod_fric(strmatch('dyhat',M_.endo_names,'exact'),:),'b-.');   
hold on
hh2=plot(timeline,100*emp_equity,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Equity payout');

%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff=figure('Name','Figure 4: Response to financial shocks only');

subplot(2,2,1)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('yhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_GDP,'g');
hold off
axis([-inf inf -14 8]);
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
ll=legend([hh1 hh2],'Model','GDP');
set(ll,'Location','NorthWest');
title('GDP');

subplot(2,2,2)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('nhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_hours,'g');
hold off
axis([-inf inf -14 8]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Hours');

subplot(2,2,3)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('byhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_debt,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Debt repurchase');

subplot(2,2,4)
hh1=plot(timeline,100*y_IRF_fin_fric(strmatch('dyhat',M_.endo_names,'exact'),:),'b-.'); 
hold on
hh2=plot(timeline,100*emp_equity,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Equity payout');
        
equity_payout_filtered_model=bpf(100*y_IRF_both_fric(strmatch('dyhat',M_.endo_names,'exact'),:)',6,32,12);
fprintf('Std(Equity payout Model): %4.3f \n',nanstd(equity_payout_filtered_model))
equity_payout_filtered_data=bpf(100*emp_equity,6,32,12);
fprintf('Std(Equity payout Data): %4.3f \n',nanstd(equity_payout_filtered_data))


%%%%%%%%%%%%%%%%%%%%%%%% Create Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff=figure('Name','Figure 5: Response to both shocks');

subplot(2,2,1)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('yhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_GDP,'g');
hold off
ll=legend('Model','GDP');
set(ll,'Location','NorthWest');
axis([-inf inf -14 8]);
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('GDP');

subplot(2,2,2)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('nhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_hours,'g');
hold off
axis([-inf inf -14 8]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Hours');

subplot(2,2,3)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('byhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_debt,'g');
hold off
axis([-inf inf -12 15]);        
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h
title('Debt repurchase');

subplot(2,2,4)
hh1=plot(timeline,100*y_IRF_both_fric(strmatch('dyhat',M_.endo_names,'exact'),:),'b-.');
hold on
hh2=plot(timeline,100*emp_equity,'g');
hold off
axis([-inf inf -12 15]);        
title('Equity payout');
plot_NBER_recessions([hh1;hh2]);     %plot shaded recession dates behind plot with handle h

fprintf('Std(Y): \t %3.2f \t %3.2f \t %3.2f \t %3.2f\n',100*std(emp_GDP), 100*std(y_IRF_prod_fric(strmatch('yhat',M_.endo_names,'exact'),:)),100*std(y_IRF_fin_fric(strmatch('yhat',M_.endo_names,'exact'),:)),100*std(y_IRF_both_fric(strmatch('yhat',M_.endo_names,'exact'),:)))
fprintf('Std(N): \t %3.2f \t %3.2f \t %3.2f \t %3.2f\n',100*std(emp_hours), 100*std(y_IRF_prod_fric(strmatch('nhat',M_.endo_names,'exact'),:)),100*std(y_IRF_fin_fric(strmatch('nhat',M_.endo_names,'exact'),:)),100*std(y_IRF_both_fric(strmatch('nhat',M_.endo_names,'exact'),:)))