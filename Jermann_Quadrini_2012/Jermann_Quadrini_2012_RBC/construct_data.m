function construct_data()
% function construct_data()
% This file constructs the data used in Jermann_Quadrini_2012_RBC.mod from the 
% raw data of JQ 2012. The correctness of the latter has been verified
% 
% Notes:
%   - running this file requires the csminwel-optimizer provided by Chris
%       Sims on his homepage (http://sims.princeton.edu/yftp/optimize/)
%   - Table 1 can only be approximately replicated as e.g. information on
%       the exact sample, the filter order, and the treatment of initial 
%       and terminal artifacts arising in the bandpass filter is missing

% Copyright (C) 2015-2016 Johannes Pfeifer
% 
%  This is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  It is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  For a copy of the GNU General Public License,
%  see <http://www.gnu.org/licenses/>.

replication=0; %set to 1 for JQ (2012) data, to 0 for Pfeifer (2016) corrected TFP

start_date=1983.75;
theta=0.36;

data_matrix = xlsread('JQ2012.xlsx',1,'I5:K238');
timeline = data_matrix(:,1);
equity_payout = data_matrix(:,2);
debt_repurchases = data_matrix(:,3);

debt_repurchases_detrended=detrend(debt_repurchases(timeline>start_date,:));
equity_payout_detrended=detrend(equity_payout(timeline>start_date,:));

total_GDP = xlsread('JQ2012.xlsx',2,'B4:B237');
survey_data = xlsread('JQ2012.xlsx',2,'G132:G237');

%% Figure 1
figure ('Name','Figure 1: Financial Flows')
h=plot(timeline,[debt_repurchases equity_payout]*100);
legend('Debt repurchase','Equity payout','Location','northwest')
axis([-inf inf -16 16]);
plot_NBER_recessions(h);     %plot shaded recession dates behind plot with handle h

%% Table 1

fprintf('Table 1\n')

equity_payout_filtered = bpf(equity_payout,6,32,12);
debt_repurchase_filtered = bpf(debt_repurchases,6,32,12);
total_GDP_filtered = bpf(total_GDP,6,32,12);

fprintf('Estimated Standard deviation of Equity: %4.3f\n', nanstd(equity_payout_filtered(timeline>start_date,:)*100));
fprintf('Estimated Standard deviation of Debt: %4.3f\n', nanstd(debt_repurchase_filtered(timeline>start_date,:)*100));
fprintf('Estimated Correlation between GDP and Equity: %4.3f\n', nancorr(equity_payout_filtered(timeline>start_date,:),total_GDP_filtered(timeline>start_date,:)));
fprintf('Estimated Correlation between GDP and Debt: %4.3f\n', nancorr(debt_repurchase_filtered(timeline>start_date,:),total_GDP_filtered(timeline>start_date,:)));

%% Construct TFP and financial shocks as described in JQ Appendix (partially based on JQ's replication files)

[numbers,text] = xlsread('JQ2012.xlsx','3-DataSet2','A2:O239');

%consistency check
if size(numbers,2)~=size(text,2)
    error('Dimensions not consistent')
end
header_rows=3;
if numbers(header_rows+1,1)~=1952
    error('Wrong number of header rows')
end
   
timeline=numbers(header_rows+1:end,1);

Nominal_Capital_Expenditures=numbers(header_rows+1:end,strcmp('Capital Expenditures',text(3,:)));
Nominal_Capital_Consumption_corporations=numbers(header_rows+1:end,strcmp('CapCon-Corp',text(3,:)));
Nominal_Capital_Consumption_noncorporations=numbers(header_rows+1:end,strcmp('CapCon-NonCo',text(3,:)));
Nominal_Net_Borrowing=numbers(header_rows+1:end,strcmp('Net Borrowing',text(3,:)));
Nominal_Business_value_added=numbers(header_rows+1:end,strcmp('Value Added',text(4,:)));
Business_Price_index=numbers(header_rows+1:end,strcmp('Price Index',text(4,:)));
Total_Private_hours=numbers(header_rows+1:end,strcmp('Hours',text(3,:)));
Real_GDP=numbers(header_rows+1:end,strcmp('Real GDP',text(3,:)));

T_periods=size(timeline,1);
Real_Capital_stock_beginning_period=NaN(T_periods+1,1);
Nominal_Debt_stock_beginning_period=NaN(T_periods+1,1);

%set starting value for perpetual inventory method
Real_Capital_stock_beginning_period(1)=21.28;
Nominal_Debt_stock_beginning_period(1)=94.12;

for t=1:T_periods
    Real_Capital_stock_beginning_period(t+1)=Real_Capital_stock_beginning_period(t)+(Nominal_Capital_Expenditures(t)-Nominal_Capital_Consumption_corporations(t)-Nominal_Capital_Consumption_noncorporations(t))*0.00025/Business_Price_index(t);
    Nominal_Debt_stock_beginning_period(t+1)=Nominal_Debt_stock_beginning_period(t)+Nominal_Net_Borrowing(t)*0.00025;
end;
Real_Capital_stock_end_period=Real_Capital_stock_beginning_period(2:end);
Nominal_Debt_stock_end_period=Nominal_Debt_stock_beginning_period(2:end);
Nominal_Capital_stock_end_period=Real_Capital_stock_end_period.*Business_Price_index;

Real_Debt_stock_end_period=Nominal_Debt_stock_end_period./Business_Price_index;
Real_Business_value_added=Nominal_Business_value_added./(Business_Price_index.*4);

capital_GDP_ratio=Real_Capital_stock_end_period./Real_Business_value_added;
debt_GDP_ratio=Real_Debt_stock_end_period./Real_Business_value_added;

log_TFP_total_GDP_end_of_period_capital=log(Real_GDP)-(1-theta)*log(Total_Private_hours)-theta*log(Real_Capital_stock_end_period);
log_TFP_total_GDP_beginning_of_period_capital=log(Real_GDP)-(1-theta)*log(Total_Private_hours)-theta*log(Real_Capital_stock_beginning_period(1:end-1));
log_TFP_total_GDP_end_of_period_capital_detrended=detrend(log_TFP_total_GDP_end_of_period_capital(timeline>start_date,:));
log_TFP_total_GDP_beginning_of_period_capital_detrended=detrend(log_TFP_total_GDP_beginning_of_period_capital(timeline>start_date,:));
log_Real_Investment_detrended=detrend(log(Nominal_Capital_Expenditures(timeline>start_date,:)./Business_Price_index(timeline>start_date,:)));

Real_Personal_Consumption=xlsread('JQ2012.xlsx','4-Data For Estimation','C4:C257');
timeline_consumption=1947:0.25:2010.25;
log_Real_Personal_Consumption=detrend(log(Real_Personal_Consumption(timeline_consumption>start_date,:)));

log_TFP_business_beginning_of_period_capital=log(Real_Business_value_added)-(1-theta)*log(Total_Private_hours)-theta*log(Real_Capital_stock_beginning_period(1:end-1));
log_TFP_business_beginning_of_period_capital_detrended=detrend(log_TFP_business_beginning_of_period_capital(timeline>start_date,:));

%% Compare various TFP measures
figure('Name','TFP Comparison')
hh=plot(timeline(timeline>start_date,:),log_TFP_total_GDP_end_of_period_capital_detrended,'b-',timeline(timeline>start_date,:),log_TFP_total_GDP_beginning_of_period_capital_detrended,'r-.',timeline(timeline>start_date,:),log_TFP_business_beginning_of_period_capital_detrended,'g--','LineWidth',1.5);
ll=legend('JQ (Real GDP, end of period capital)','Real GDP, beginning of period capital','Business Value Added, beginning of period capital');
set(ll,'Location','NorthWest')
axis tight
plot_NBER_recessions(hh);     %plot shaded recession dates behind plot with handle h

log_Real_Business_value_added_detrended=detrend(log(Real_Business_value_added(timeline>start_date,:)));
log_Real_GDP_detrended=detrend(log(Real_GDP(timeline>start_date,:)));
log_Total_Private_hours_detrended=detrend(log(Total_Private_hours(timeline>start_date,:)));

log_Real_Capital_stock_end_period_detrended=detrend(log(Real_Capital_stock_end_period(timeline>start_date,:)));
log_Real_Debt_stock_end_period_detrended=detrend(log(Real_Debt_stock_end_period(timeline>start_date,:)));
%% construct financial friction xi
xi = -1.5489*log_Real_Capital_stock_end_period_detrended+0.5489*log_Real_Debt_stock_end_period_detrended+log_Real_Business_value_added_detrended;

%% Estimate VAR(1)
if replication
    x_orig = [log_TFP_total_GDP_end_of_period_capital_detrended xi]';
    x= x_orig;
else
    x = [log_TFP_business_beginning_of_period_capital_detrended xi]';
    x_orig = [log_TFP_total_GDP_end_of_period_capital_detrended xi]';
    Ymat_orig = x_orig(:,2:end)';
    Zmat_orig = x_orig(:,1:end-1)';

    Z_orig=Zmat_orig';
    resid_orig=(eye(size(Z_orig,2))-Z_orig'/(Z_orig*Z_orig')*Z_orig)*Ymat_orig;
end

Ymat = x(:,2:end)';
Zmat = x(:,1:end-1)';
timeline_resids=timeline(timeline>start_date);
timeline_resids=timeline_resids(2:end);
[nobs,nvars]=size(Ymat);
Z=Zmat';
ymat=Ymat';
bhat=kron((Z*Z')\Z,eye(nvars))*ymat(:);
display('The estimate of the A matrix is')
bhat_OLS=reshape(bhat,nvars,nvars)          

%% Plot shocks
resid=(eye(size(Z,2))-Z'/(Z*Z')*Z)*ymat';

figure
if ~replication
    h1=scatter(resid_orig(:,1),resid_orig(:,2),'+','LineWidth',0.5);
    hold on
    h2=scatter(resid(:,1),resid(:,2));
    xlabel('$\varepsilon_{z,t}$','Interpreter','latex','FontSize',18)
    ylabel('$\varepsilon_{\xi,t}$','Interpreter','latex','FontSize',18)
    hh=legend([h1;h2],'JQ (Real GDP, end of period capital)','Business Value Added, beginning of period capital');
else
    h1=scatter(resid(:,1),resid(:,2));
    xlabel('$\varepsilon_{z,t}$','Interpreter','latex','FontSize',18)
    ylabel('$\varepsilon_{\xi,t}$','Interpreter','latex','FontSize',18)
    hh=legend(h1,'JQ (Real GDP, end of period capital)');
end
set(hh,'Location','NorthOutside')

figure('Name','Time series of shocks')
h=plot(timeline_resids,resid(:,1),timeline_resids,resid(:,2));
plot_NBER_recessions(h);
hh=legend(h,'$\varepsilon_{z,t}$','$\varepsilon_{\xi,t}$');
set(hh,'Interpreter','latex','FontSize',18,'Location','SouthWest')
axis tight

disp(['Estimated sigma_z is ' num2str(std(resid(:,1)))]); 
disp(['Estimated sigma_xi is ' num2str(std(resid(:,2)))]);
cov_matrix=cov(resid);
[rho,pval]=corr(resid(:,1),resid(:,2))
disp(['Estimated covariance between the two shocks is ' num2str(cov_matrix(2,1))]);
disp(['Estimated correlation between the two shocks is ' num2str(cov_matrix(2,1)/((std(resid(:,1)))*(std(resid(:,2)))))])
if replication
    save('innovations_replication.mat','resid','bhat_OLS','cov_matrix')         % save residuals for .mod-file
else    
    save('innovations_corrected.mat','resid','bhat_OLS','cov_matrix')         % save residuals for .mod-file    
end

%% test restriction that covariance matrix is diagonal

L_OLS=chol(cov_matrix,'lower');

AR_pars=nvars^nvars;
total_pars=AR_pars+nvars*(nvars+1)/2;

x_start=zeros(total_pars,1);
x_start(1:AR_pars)=bhat_OLS(:);
x_start(AR_pars+1:total_pars)=cholmat2vec(L_OLS);

H0=1e-2*eye(total_pars);
crit=1e-10;
nit=1000;

par_indices_to_optimize=(1:total_pars);
%% run ML estimation
[fhat,xhat]=csminwel(@VAR_ML_restricted,x_start,H0,[],crit,nit,x_start,par_indices_to_optimize,ymat,Z);

%% run restricted ML estimation
L_OLS=chol(diag(diag(cov_matrix)),'lower');
x_restricted_fixed_par=xhat;
x_restricted_fixed_par(AR_pars+1:total_pars)=cholmat2vec(L_OLS);
par_indices_to_optimize=(1:total_pars);
par_indices_to_optimize(AR_pars+2)=[];
x_start=x_restricted_fixed_par;
x_start(AR_pars+2)=[];
H0=1e-4*eye(length(x_start));
[fhat_restricted,xhat_restricted]=csminwel(@VAR_ML_restricted,x_start,H0,[],crit,nit,x_restricted_fixed_par,par_indices_to_optimize,ymat,Z);

fprintf('Likelihood ratio test: %4.3f\n',1-chi2cdf(-2*(fhat-fhat_restricted),1))
%One can use Wolfram Alpha to evaluate CDF[chisquared[1],74]

%reconstruct parameters
betta=reshape(xhat(1:AR_pars,1),nvars,...
    nvars);
L=vec2cholmat(xhat(AR_pars+1:end,1));
Sigma_U=L*L';

%reconstruct parameters
x_restricted_full_vector=x_restricted_fixed_par;
x_restricted_full_vector(par_indices_to_optimize,:)=xhat_restricted;

betta_restricted=reshape(x_restricted_full_vector(1:AR_pars,1),nvars,...
    nvars);
L_restricted=vec2cholmat(x_restricted_full_vector(AR_pars+1:end,1));
Sigma_U_restricted=L_restricted*L_restricted';


figure ('Name','Figure 2: Time Series of Shocks to Productivity and Financial Conditions');
subplot(2,3,1);
plot(timeline(timeline>1983.75),x(1,:)*100);
axis([-inf inf -8 8]);
title('Level of productivity, z');

subplot(2,3,2);
plot(timeline(timeline>1983.75),xi*100);
axis([-inf inf -8 8]);
title('Level of fin. conditions, xi');

subplot(2,3,4);
plot(timeline(timeline>1984),resid(:,1)*100);
axis([-inf inf -3 3]);
title('Innovations to prod.');

subplot(2,3,5);
plot(timeline(timeline>1984),resid(:,2)*100);
axis([-inf inf -3 3]);
title('Innovations to fin. cond.');

subplot(2,3,6);
plot(timeline(timeline>1984),-resid(:,2)*100);
axis([-inf inf -3 4]);
hold on;
plot(timeline(timeline>1984),[NaN(23,1);survey_data/26],'--r');        % need to rescale the survey data to make it fit the paper version
hold off;
title('Index of tightening stand.');

%% Prepare detrended data for counterfactual exercises
emp_data_JQ_original=xlsread('JQ2012.xlsx',2,'B133:F237');

log_Real_Investment_detrended=log_Real_Investment_detrended(2:end);
log_Real_Business_value_added_detrended=log_Real_Business_value_added_detrended(2:end);
log_Real_Personal_Consumption=log_Real_Personal_Consumption(2:end);
log_Real_GDP_detrended=log_Real_GDP_detrended(2:end);
log_Total_Private_hours_detrended=log_Total_Private_hours_detrended(2:end);
equity_payout_detrended=equity_payout_detrended(2:end);
debt_repurchases_detrended=debt_repurchases_detrended(2:end);

data_timeline=timeline(timeline>1984);

save('emp_data.mat','data_timeline','log_Real_Business_value_added_detrended','log_Real_Investment_detrended',...
    'log_Real_Personal_Consumption','log_Real_GDP_detrended','debt_repurchases_detrended','equity_payout_detrended',...
    'log_Total_Private_hours_detrended');