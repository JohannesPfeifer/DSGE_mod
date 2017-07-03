function [IRF_mat,IRF_weighting,IRF_quantiles]=get_empirical_IRFs(IRF_length,parametric_dummy)
% function [IRF_mat,IRF_weighting,IRF_quantiles]=get_empirical_IRFs(IRF_length)
% Estimates the IRFs to a goverment spending shock following
% Blanchard/Perotti (2002), employing a residual bootstrap to derive the
% confidence bands; as in Blanchard/Perotti, the confidence bands are
% symmetric around the mean response, assuming a normal distribution. Alternatively, the 
% respective quantiles of the distribution can be used. In principle, one
% may additionally use a bias correction (e.g. Kilian (1998)).
% Unstable draws are discarded in the bootstrapping.
%
% Inputs:
%   IRF_length          [scalar]                IRF horizon
%   IRF_target          [nperiods by nvars]     matrix of target IRFs
%   parametric_dummy    [scalar]                1: use symmetric bands based on normal distribution
%                                               0: use quantiles
% Outputs:
%   IRF_mat             [nperiods by nvars]  matrix of empirical IRFs
%   IRF_weighting       [nperiods*nvars by nperiods*nvars]     matrix of weights for IRF matching 
%                                                               (a matrix with the inverse of the 
%                                                               variances of the IRFS on the diagonal)
%   IRF_quantiles       [nvars by nperiods by 2]  IRF quantiles from bootstrapping
%
% Copyright (C) 2016-17 Benjamin Born and Johannes Pfeifer
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.

if nargin<2
    parametric_dummy=1;
end

%use 4 lags, constant, and deterministic trend
number_of_lags = 4;
constant_dummy =1;
trend_dummy = 1;

%% Data preparation
real_gdp = xlsread('BP.xlsx',1,'B8:B285');
gdp_deflator = xlsread('BP.xlsx',1,'D8:D285');
nom_gov_cons_inv = xlsread('BP.xlsx',1,'F8:F285');
nom_priv_cons_nd = xlsread('BP.xlsx',1,'J8:J285');
nom_priv_cons_serv = xlsread('BP.xlsx',1,'L8:L285');
nipa_pop=xlsread('BP.xlsx','N8:N285');

% construct real per capita values
Y = real_gdp./nipa_pop;
G = nom_gov_cons_inv./nipa_pop./gdp_deflator;
C = (nom_priv_cons_nd+nom_priv_cons_serv)./nipa_pop./gdp_deflator;

timeline = (1947:0.25:2016.25)';

gov_spend_to_gdp_share=nom_gov_cons_inv./(real_gdp.*gdp_deflator/100);
gov_spend_to_gdp_share_mean=mean(gov_spend_to_gdp_share(timeline>1953.75 & timeline<2008,:));

cons_to_gdp_share=(nom_priv_cons_nd+nom_priv_cons_serv)./(real_gdp.*gdp_deflator/100);
cons_to_gdp_share_mean=mean(cons_to_gdp_share(timeline>1953.75 & timeline<2008,:));


x =[log(G) log(Y) log(C)];
var_names={'G','Y','C'};
G_pos=1;
Y_pos=2;
C_pos=3;

%restrict to sample before financial crisis
x=x(timeline>1953.75 & timeline<2008,:);

%create regressor matrix of lagged variables
Z=[];

for lagnumber=1:number_of_lags
    tempmat = lagmatrix(x,lagnumber);
    Z = [Z tempmat(number_of_lags+1:end,:)];
end
[T, nvars] =size(x);

%add trend and constant
if trend_dummy==1
    Z = [(1:T-number_of_lags)' Z];
end

if constant_dummy==1
    Z = [ones(T-number_of_lags,1) Z];
end

%get OLS estimates using the Kronecker formula (e.g. Luetkepohl (2005))
Z = Z';
Y=x(number_of_lags+1:end,:)';
T=size(Z,2);
bhat = kron((Z*Z')\Z,eye(nvars))*Y(:);
betta_OLS=reshape(bhat,nvars,nvars*number_of_lags+constant_dummy+trend_dummy);

% get residuals and compute covariance matrix
resids=(eye(T)-Z'/(Z*Z')*Z)*Y';
L_OLS=chol(cov(resids),'lower');

%construct IRFs using the companion form, i.e. transforming the model into a
%VAR(1)
[F,Q]=build_companionform(betta_OLS(:,constant_dummy+trend_dummy+1:end),number_of_lags,...
    L_OLS);

n_augmented_vars=size(F,1);
IRFs_point=zeros(n_augmented_vars,IRF_length); %initialize IRF matrix

%set impulse to 1 percent of GDP
shock_vector=zeros(n_augmented_vars,1);
shock_pos=G_pos;
shock_vector(shock_pos)=1/Q(shock_pos,shock_pos); %set to 1 percent of G

% generate IRFs
for boot_iter=1:IRF_length
    IRFs_point(:,boot_iter)=F^(boot_iter-1)*Q*shock_vector;
end

%select IRFs for matching
IRF_mat=[IRFs_point(G_pos,:)' IRFs_point(Y_pos,:)'];

%% get uncertainty bands
nbootstraps = 100;

generatedseries=NaN(T+number_of_lags,nvars); %initialize matrix
IRFs=zeros(n_augmented_vars,IRF_length,nbootstraps);

generatedseries(1:number_of_lags,:)=x(1:number_of_lags,:); %start with true starting values
determmatrix = Z(1:constant_dummy+trend_dummy,:); %select constant and trend if necessary

for boot_iter=1:nbootstraps
    unstable_indicator=1;
    iter=1;
    while unstable_indicator==1 && iter<1000
        iter=iter+1;
        indices=ceil(rand(T,1)*T); %randomly generate index series to pick out residuals; scale the [0,1]-distribution up with nobs
        bootresids=resids(indices,:); %select residuals

        for jj=number_of_lags+1:T+number_of_lags %recursively generate new time series
            ylagvectemp=generatedseries(jj-number_of_lags:jj-1,:); %stack endogenous variables in format compatible with matrix multiplication
            ylagvectemp=ylagvectemp(number_of_lags:-1:1,:); %reverse order for vectorization
            ylagvector=vec(ylagvectemp');     
            generatedseries(jj,:) = determmatrix(:,jj-number_of_lags)'*betta_OLS(:,1:constant_dummy+trend_dummy)' + (betta_OLS(:,1+constant_dummy+trend_dummy:end)*ylagvector)'+ bootresids(jj-number_of_lags,:);
        end
        % construct lagged regressor matrix
        Z_generated=[];
        for lagnumber=1:number_of_lags
            tempmat = lagmatrix(generatedseries,lagnumber);
            Z_generated = [Z_generated tempmat(number_of_lags+1:end,:)];
        end
        Z_reg_mat = [determmatrix; Z_generated']; %add deterministics
        Y_generated=generatedseries(number_of_lags+1:end,:)'; %generate dependent variables
        % reestimate coefficients
        bhat_boot = kron((Z_reg_mat*Z_reg_mat')\Z_reg_mat,eye(nvars))*Y_generated(:);
        betta_boot=reshape(bhat_boot,nvars,nvars*number_of_lags+constant_dummy+trend_dummy);
        resids_boot=(eye(T)-Z_reg_mat'/(Z_reg_mat*Z_reg_mat')*Z_reg_mat)*Y_generated';
        L_boot=chol(cov(resids_boot),'lower');

        %generate IRFs
        [F_boot,Q_boot]=build_companionform(betta_boot(:,constant_dummy+trend_dummy+1:end),number_of_lags,L_boot);
        if max(abs(eig(F_boot)))>0.999 %discard unstable draws
            unstable_indicator=1;
            continue
        else
            unstable_indicator=0;
        end

        shock_vector=zeros(size(F_boot,1),1);
        shock_vector(shock_pos)=1/Q_boot(shock_pos,shock_pos); %set to 1 percent of government spending

        for irf_iter=1:IRF_length
            IRFs(:,irf_iter,boot_iter)=F_boot^(irf_iter-1)*Q_boot*shock_vector;
        end
    end
end
if iter==1000
    error('Was not able to generate a stable draw in 1000 tries')
end

%get pointwise variance of IRFs across draws and compute weighting matrix
IRF_variances=var(IRFs([1,3],:,:),0,3)';
IRF_weighting=inv(diag(IRF_variances(2:end)));

% compute confidence bands for plots
if parametric_dummy
    IRF_std=std(IRFs,0,3);
    IRF_quantiles(:,:,1)=IRFs_point-1.96*IRF_std;
    IRF_quantiles(:,:,2)=IRFs_point+1.96*IRF_std;
else
    IRF_quantiles = quantile(IRFs(1:nvars,:,:),[0.025 0.975],3);    
end

% plot results, measured in percent of GDP
figure
subplot(1,3,1);
plot(1:IRF_length,IRF_quantiles(1,:,1),'r--',1:IRF_length,IRF_quantiles(1,:,2),'r--',1:IRF_length,IRFs_point(1,:),'b-')
title(var_names(G_pos));
ylim([-1 2.5])
subplot(1,3,2);
plot(1:IRF_length,1/gov_spend_to_gdp_share_mean*IRF_quantiles(2,:,1),'r--',1:IRF_length,1/gov_spend_to_gdp_share_mean*IRF_quantiles(2,:,2),'r--',1:IRF_length,1/gov_spend_to_gdp_share_mean*IRFs_point(2,:),'b-')
title(var_names(Y_pos));
ylim([-1 2.5])

subplot(1,3,3);
plot(1:IRF_length,cons_to_gdp_share_mean/gov_spend_to_gdp_share_mean*IRF_quantiles(3,:,1),'r--',1:IRF_length,cons_to_gdp_share_mean/gov_spend_to_gdp_share_mean*IRF_quantiles(3,:,2),'r--',1:IRF_length,cons_to_gdp_share_mean/gov_spend_to_gdp_share_mean*IRFs_point(3,:),'b-')
title(var_names(C_pos));
ylim([-1 2.5])

if nargout>2
    IRF_quantiles=[IRF_quantiles([G_pos, Y_pos],:,:)];
end

