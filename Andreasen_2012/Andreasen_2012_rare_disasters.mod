/*
 * This file replicates Table 2 (Benchmark, Non-symmetric Case I and Case II)
 * and Table 4 (Benchmark) of the NK model with Epstein-Zin preferences and
 * asymmetric shocks as in:
 * Andreasen (2012): "On the effects of rare disasters and uncertainty shocks for risk premia in non-linear DSGE models",
 * Review of Economic Dynamics, 15, pp. 295-316.
 *
 * Notes:
 *  - Dynare's perturbation solver assumes Gaussian innovations. Therefore,
 *    in order to do simulations with non-symmetric shocks one has to manually
 *    re-compute some perturbation matrices as the third-order product moments
 *    are assumed to be zero and hence oo_.dr.ghs3 is set to zero automatically.
 *  - If the second-moments do not correspond to what is declared in the
 *    shocks block (i.e. M_.Sigma_e), then one also needs to re-compute
 *    oo_.dr.ghs2, oo_.dr.ghxss, oo_.dr.ghuss.
 *  - The script perturbation_solver_nonsymmetric_order3.m re-computes
 *    oo_.dr.ghs3, and possibly oo_.dr.ghs2, oo_.dr.ghxss, oo_.dr.ghuss,
 *    provided the vectors of second-order product moments SIGMA2 = E[kron(u,u)]
 *    and third-order product moments SIGMA3 = E[kron(u,kron(u,u))]
  * - The model is based on the replication codes of Martin Andreasen.
 *    Note that the original implementation uses the SGU model framework
 *    and a perturbation solver tailored specific to the model.
 *  - This mod file was created with Dynare 5.2
 *  - The runtime is a couple of minutes depending on the number of periods 
 *
 * Additional functions:
 *  - perturbation_solver_nonsymmetric_order3.m
 *    re-computes the perturbation matrices in case of non-Gaussian and
 *    assymetric shocks up to a third-order approximation
 *  - simult_nonsymmetric_order3.m: does stochastic simulations given the
 *    third-order (unpruned) state-space system taking ghs3 into account
 *
 * This implementation was written by Willi Mutschler (@wmutschl).
 * If you spot mistakes, email me at willi@mutschler.eu
 *
 * Please note that the following copyright notice only applies to this Dynare
 * implementation of the model.
 */

/*
 * Copyright © 2022 Willi Mutschler
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
 * see <https://www.gnu.org/licenses/>.
 */

/*
 **************************************************************************
 * Declarations based on:
 * - Run_Rotemberg_Model_growth.m
 * - NK_Rotemberg_model.m
 * - Anal_SDF_Rotemberg.m
 * - Anal_SDF_neutral.m
 * - Get_Bond_Prices_Mom3_Log.m
 * - Get_Data_sim.m
 **************************************************************************
 */

var
  ln_r                       (long_name='log of short-term nominal interest rate')
  ln_g                       (long_name='log of government spending')
  ln_p1                      (long_name='log of short-term bond price')
  ln_y                       (long_name='log of output')
  ln_c                       (long_name='log of consumption')
  ln_n                       (long_name='log of labor supply')
  ln_evf                     (long_name='log of expected value of the value function')
  ln_pai                     (long_name='log of gross inflation')
  @#for j in 2:40
  ln_p@{j}                   (long_name='log of bond price for maturity @{j} under Rotemberg pricing kernel')
  @#endfor
  @#for j in 2:40
  ln_q@{j}                   (long_name='log of bond price for maturity @{j} under neutral pricing kernel')
  @#endfor
  @#for j in 2:40
  ln_tp@{j}                  (long_name='log of term premium for maturity @{j}')
  @#endfor

  //reporting variables (we will declare them as VAROBS below for simpler access)
  Gr_C       ${\Delta c}$    (long_name='Δc_t: growth rate in consumption: annualized and in pct')
  Infl       ${\pi}$         (long_name='π_t: inflation rate: annualized and in pct')
  R1         ${r}$           (long_name='r_t: short-term interest rate: annualized and in pct')
  R40        ${r_{40}}$      (long_name='r40_t: long-term interest rate: annualized and in pct')
  Slope      ${r_{40}-r}$    (long_name='r40_t-r_t: Slope of yield curve')
  TP         ${P_{40}}$      (long_name='P40_t: term premia: annualized basis points')
  xhr40      ${xhr_{40}}$    (long_name='xhr40_t: 10 year excess holding period return')
  varsdf     ${Var(M_{t,t+1})}$ (long_name='variance of the nominal stochastic discount factor')
  ln_a       ${a_t}$         (long_name='log of technology level')
;

varobs Gr_C Infl R1 R40 Slope TP xhr40 varsdf ln_a;

varexo
  epsA       ${\epsilon_a}$    (long_name='productivity shock (might be non-Gaussian)')
  epsG       ${\epsilon_g}$    (long_name='government spending shock (Gaussian)')
  epsR       ${\epsilon_r}$    (long_name='monetary polilcy shock (Gaussian)')
;

parameters
  GAMA       ${\gamma}$           (long_name='first parameter that determines intertemporal elasticity of substitution')
  NU         ${\nu}$              (long_name='second parameter that determines intertemporal elasticity of substitution')
  BETTA      ${\beta}$            (long_name='discount factor')
  ALFA       ${\alpha}$           (long_name='degree of relative risk-aversion')
  THETA      ${\theta}$           (long_name='capital elasticity in Cobb-Douglas production function')
  ETA        ${\eta}$             (long_name='love of variety parameter in Dixit-Stiglitz aggregator')
  XI         ${\xi}$              (long_name='Rotemberg quadratic price adjustment cost coefficient')
  DELTA      ${\delta}$           (long_name='depreciation rate of capital')
  PHI_PAI    ${\phi_\pi}$         (long_name='inflation feedback Taylor Rule')
  PHI_Y      ${\phi_y}$           (long_name='output gap feedback Taylor Rule')
  RHOA       ${\rho_a}$           (long_name='persistence parameter productivity process')
  RHOG       ${\rho_g}$           (long_name='persistence parameter government spending process')
  RHOR       ${\rho_r}$           (long_name='persistence parameter Taylor rule')
  N          ${n_{ss}}$           (long_name='steady-state labor supply')
  G_O_Y      ${g_{ss}/y_{ss}}$    (long_name='steady-state ratio of government spending to output')
  PAI        ${\pi_{ss}}$         (long_name='inflation target')
  MUZ        ${\mu_{z,ss}}$       (long_name='deterministic trend')
  Kss        ${\bar{k}}$          (long_name='steady-state capital stock')
  AA                              (long_name='normalization constant equal to minus the steady-state of expected value function')
  Ass                             (long_name='steady-state level of technology')
;

/*
 **************************************************************************
 * Calibration based on
 * - Table 1
 * - Run_Rotemberg_Model_growth.m
 * - NK_Rotemberg_model_ss.m
 **************************************************************************
 */

Ass = 1;
GAMA = 2.5;
NU = 0.35;
BETTA = 0.9995;
ALFA = -110;
THETA = 0.36;
ETA = 6;
ALFAp = 0.75;
DELTA = 0.025;
PHI_PAI = 1.5;
PHI_Y = 0.30;
RHOA = 0.98;
RHOG = 0.90;
RHOR = 0.85;
N = 0.38;
G_O_Y = 0.17;
PAI = 1.008;
MUZ = 1.005;
XI = (1-THETA+ETA*THETA)*(ETA-1)*ALFAp/((1-ALFAp)*(1-THETA)*(1-ALFAp*BETTA*MUZ^(NU*(1-GAMA)))); % in Table 1 it is given by 260; however, in the replication codes it is an endogenous parameter
IES = 1/(1-NU*(1-GAMA));
RRA = (GAMA+ALFA*(1-GAMA));
fprintf('Chosen calibration implies:\n');
fprintf('  - Inverse Elasticity of Substitution (IES): %f\n',IES);
fprintf('  - Risk Aversion (RRA): %f\n',RRA);

/*
 **************************************************************************
 * Model equations based on
 * - Run_Rotemberg_Model_growth.m
 * - NK_Rotemberg_model.m
 * - Anal_SDF_Rotemberg.m
 * - Anal_SDF_neutral.m
 * - Get_Bond_Prices_Mom3_Log.m
 * - Get_Data_sim.m
 **************************************************************************
 */

model;
    //do exp transform on most variables (not all) such that model equations contain actual variables and not logged variables
    #r_back = exp(ln_r(-1));
    #a_back = exp(ln_a(-1));
    #g_back = exp(ln_g(-1));

    #a_curr = exp(ln_a);
    #g_curr = exp(ln_g);
    #p1_curr = exp(ln_p1);
    #r_curr = exp(ln_r);
    #y_curr = exp(ln_y);
    #c_curr = exp(ln_c);
    #n_curr = exp(ln_n);
    #evf_curr = exp(ln_evf);
    #pai_curr = exp(ln_pai);

    #a_fwrd = exp(ln_a(+1));
    #g_fwrd = exp(ln_g(+1));
    #p1_fwrd = exp(ln_p1(+1));
    #r_fwrd = exp(ln_r(+1));
    #y_fwrd = exp(ln_y(+1));
    #c_fwrd = exp(ln_c(+1));
    #n_fwrd = exp(ln_n(+1));
    #evf_fwrd = exp(ln_evf(+1));
    #pai_fwrd = exp(ln_pai(+1));

    #varsdf_curr = varsdf; //no exp transform
    #varsdf_fwrd = varsdf(+1); //no exp transform

    //auxiliary expressions and parameters
    #muz_curr = MUZ;
    #muz_fwrd = MUZ;
    #PAIss = exp(steady_state(ln_pai));
    #Rss = exp(steady_state(ln_r));
    #Yss = exp(steady_state(ln_y));
    #Gss = G_O_Y*Yss;
    #w_curr = (1-NU)*c_curr/(NU*(1-n_curr)); //FOC for household with respect to labor
    #mc_curr = w_curr/((1-THETA)*a_curr*Kss^THETA*n_curr^(-THETA)); //FOC for firms with respect to labor
    #mvf_fwrd = -(c_fwrd^(NU*(1-GAMA))*(1-n_curr)^((1-NU)*(1-GAMA))/(1-GAMA)-BETTA*muz_fwrd^(NU*(1-GAMA))*AA*evf_fwrd^(1/(1-ALFA))); //The expression for minus the value function (use the fact that the trend is deterministic)
    #mu_la_fwrd = (AA*evf_curr^(1/(1-ALFA))/mvf_fwrd)^ALFA*muz_fwrd^(NU*(1-GAMA)-1)*(c_fwrd/c_curr)^(NU*(1-GAMA)-1)*((1-n_fwrd)/(1-n_curr))^((1-NU)*(1-GAMA)); //The ratio of lambda_fwrd/lambda_curr

    //actual model equations
    [name='Expected value of the value function']
    evf_curr = (mvf_fwrd/AA)^(1-ALFA);
    [name='The one period bond price']
    p1_curr = BETTA*mu_la_fwrd*1/pai_fwrd;
    [name='The one period interest rate']
    1 = BETTA*mu_la_fwrd*r_curr/pai_fwrd;
    [name='The FOC firms with respect to prices']
    mc_curr = (ETA-1)/ETA - BETTA*(AA*evf_curr^(1/(1-ALFA))/mvf_fwrd)^ALFA*muz_fwrd^(NU*(1-GAMA))*(c_fwrd/c_curr)^(NU*(1-GAMA)-1)*((1-n_fwrd)/(1-n_curr))^((1-NU)*(1-GAMA))*XI/ETA*(pai_fwrd/PAIss-1)*pai_fwrd*y_fwrd/(PAIss*y_curr) + XI/ETA*(pai_curr/PAIss-1)*pai_curr/PAIss;
    [name='The Taylor rule']
    log(r_curr/Rss) = RHOR*log(r_back/Rss)+ PHI_PAI*log(pai_curr/PAIss) + PHI_Y*log(y_curr/Yss) + epsR;
    [name='The output level']
    y_curr = a_curr*Kss^THETA*n_curr^(1-THETA);
    [name='The budget resource constraint']
    y_curr = c_curr + g_curr + DELTA*Kss;
    [name='Technology shocks']
    log(a_curr/Ass) = RHOA*log(a_back/Ass) + epsA;
    [name='Shocks to government spendings']
    log(g_curr/Gss) = RHOG*log(g_back/Gss) + epsG;

    //compute bond prices up to a given maturity using the log-transformation
    #M_Q = 1/r_curr; //pricing kernel given risk neutral pricing (Q)
    #M_P = BETTA*mu_la_fwrd/pai_fwrd; //pricing kernel given Rotemberg prices(P)
    #q1_curr = p1_curr;
    #q1_fwrd = p1_fwrd;
    @#for j in 2:40
    #p@{j}_curr = exp(ln_p@{j});
    #p@{j}_fwrd = exp(ln_p@{j}(+1));
    #q@{j}_curr = exp(ln_q@{j});
    #q@{j}_fwrd = exp(ln_q@{j}(+1));
    [name='The yield curve: p@{j}']
    p@{j}_curr = M_P*p@{j-1}_fwrd;
    [name='The yield curve: q@{j}']
    q@{j}_curr = M_Q*q@{j-1}_fwrd;
    [name='The term premium: ln_tp@{j}']
    ln_tp@{j} = -1/@{j}*(ln_p@{j}-ln_q@{j});
    @#endfor

    //auxiliary reporting equations
    [name='Δc_t: growth rate in consumption: annualized and in pct']
    Gr_C = 400*(log(MUZ) + ln_c - ln_c(-1));
    [name='π_t: inflation rate: annualized and in pct']
    Infl = 400*ln_pai;
    [name='r_t: short-term interest rate: annualized and in pct']
    R1 = 400*ln_r;
    [name='r40_t: long-term interest rate: annualized and in pct']
    R40 = 400*(-1/40*ln_p40);
    [name='r40_t-r_t: Slope of yield curve']
    Slope = R40 - R1;
    [name='P40_t: term premia: annualized basis points']
    TP = 400*ln_tp40;
    [name='xhr40_t: 10 year excess holding period return']
    xhr40 = (ln_p39 - ln_p40(-1) - ln_r(-1))*400;
    [name='The variance of the nominal stochastic discount factor']
    varsdf_curr = M_P^2 - M_Q^2;
end;


/*
 **************************************************************************
 * Steady State based on
 * - NK_Rotemberg_model_ss.m
 * - Anal_SDF_Rotemberg.m
 * - Anal_SDF_neutral.m
 **************************************************************************
 */

steady_state_model;
    PSS = BETTA*MUZ^(NU*(1-GAMA)-1)/PAI;
    RSS = 1/PSS;
    MCSS = (ETA-1)/ETA;
    C_O_Y = ((1-NU)*N/(NU*(1-N)*MCSS*(1-THETA)))^-1;
    K_O_Y = (1-C_O_Y-G_O_Y)/DELTA;
    N_O_Y = (1/Ass*K_O_Y^(-THETA))^(1/(1-THETA));
    WSS = MCSS*(1-THETA)*(N_O_Y)^-1;
    CSS = WSS*NU*(1-N)/(1-NU);
    YSS = C_O_Y^-1*CSS;
    Kss = K_O_Y*YSS; %endogenous parameter (updated after call to steady)
    MVFSS = CSS^(NU*(1-GAMA))*(1-N)^((1-NU)*(1-GAMA))/((1-GAMA)*(BETTA*MUZ^(NU*(1-GAMA))-1));
    AA  = MVFSS; %endogenous parameter (updated after call to steady)
    EVFSS= (MVFSS/AA)^(1-ALFA);
    GSS = G_O_Y*YSS;
    MQ = 1/RSS;
    MP = BETTA*(AA*EVFSS^(1/(1-ALFA))/MVFSS)^ALFA*MUZ^(NU*(1-GAMA)-1)*(CSS/CSS)^(NU*(1-GAMA)-1)*((1-N)/(1-N))^((1-NU)*(1-GAMA))/PAI;
    ln_p = log(PSS);
    ln_r = log(RSS);
    ln_y = log(YSS);
    ln_c = log(CSS);
    ln_n = log(N);
    ln_evf = log(EVFSS);
    ln_pai = log(PAI);
    ln_a = log(Ass);
    ln_g = log(GSS);
    @#for j in 1:40
    ln_p@{j} = ln_p*@{j};
    ln_q@{j} = ln_p*@{j};
    ln_tp@{j} = -1/@{j}*(ln_p@{j}-ln_q@{j});
    @#endfor
    Gr_C = 400*log(MUZ);
    Infl = 400*log(PAI);
    R1 = 400*ln_r;
    R40 = 400*(-1/40*ln_p40);
    Slope = R40 - R1;
    TP = 400*ln_tp40;
    xhr40 = (ln_p39 - ln_p40 - ln_r)*400;
    varsdf = 0;
end;

/*
 **************************************************************************
 * shock variance specification based on
 * - Table 1
 * - Run_Rotemberg_Model_growth.m
 * to speed up computations set the same covariance matrix as in SIGMA2 below
 * as perturbation_solver_nonsymmetric_order3.m then does not need to re-compute
 * ghs2, ghxss and ghuss if SIGMA2 is equal to M_.Sigma_e
 **************************************************************************
 */
shocks;
    %var epsA = 1;
    %var epsG = 1;
    %var epsR = 1;    
    var epsA = 0.0075^2;
    var epsG = 0.004^2;
    var epsR = 0.003^2;
end;

/*
 **************************************************************************
 * pre-computations assuming Gaussian shocks with covariance M_.Sigma_e
 **************************************************************************
 */
steady; //important to compute steady-state first as this updates the endogenous parameters
stoch_simul(order=3,irf=0,periods=0,nocorr,nodecomposition,nofunctions,nomoments); //third-order perturbation matrices from symmetric Gaussian shocks


/*
 **************************************************************************
 * do manual simulations
 **************************************************************************
 */
options_.periods = 2e6; % value used in replication codes
fprintf('Please wait, depending on options_.periods this might take a couple of minutes\n~3 minutes in total for options_.periods = 2e6 on a MacBook Air M1 (2021) with Dynare 5.2 and MATLAB R2022b (under Rosetta)\n')

% shock variance specification based on Table 1 and Run_Rotemberg_Model_growth.m
SIGMA2 = zeros(M_.exo_nbr,M_.exo_nbr);
SIGMA2(find(ismember(M_.exo_names,'epsA')),find(ismember(M_.exo_names,'epsA'))) = 0.0075^2;
SIGMA2(find(ismember(M_.exo_names,'epsG')),find(ismember(M_.exo_names,'epsG'))) = 0.004^2;
SIGMA2(find(ismember(M_.exo_names,'epsR')),find(ismember(M_.exo_names,'epsR'))) = 0.003^2;

% draw shocks
randn('seed',1) % this is the seed used in Run_Rotemberg_Model_growth.m
exo = transpose(chol(SIGMA2)*randn(M_.exo_nbr,options_.periods)); % draw from standard normal and multiply with standard deviations

% Benchmark: all shocks are Gaussian
if isequal(SIGMA2,M_.Sigma_e)
    oo_.dr.ghs3 = zeros(M_.endo_nbr,1);
else
    oo_ = perturbation_solver_nonsymmetric_order3(M_,oo_,SIGMA2(:),SIGMA3(:));
end
RESULTS.BENCHMARK.exo = exo;
fprintf('Simulating Benchmark model\n')
RESULTS.BENCHMARK.y = simult_nonsymmetric_order3(M_,oo_,options_,exo); % simulate data for BENCHMARK model

% Non-Symmetric: shocks are drawn from a mixed-normal
for model_variant = ["CASEI" "CASEII"]
    SIGMA3 = zeros(M_.exo_nbr,M_.exo_nbr,M_.exo_nbr);
    if model_variant == "CASEI"
        PHI = -38.6667; % this is the value used in Table 2 Non-Symmetric Case I
    elseif model_variant == "CASEII"
        PHI = -7.00; % this is the value used in Table 2 Non-Symmetric Case II
    end
    randn('seed',10); rand('seed',11); % these are the seeds used in rand_mixN.m
    tmpu = rand(options_.periods,1);     tmp1 = randn(options_.periods,1);    tmp2 = repmat(PHI,options_.periods,1) + randn(options_.periods,1)*0.00001;
    epsA_skew = zeros(options_.periods,1);
    for t=1:options_.periods
       if tmpu(t,1) < 0.9958
           epsA_skew(t,1) = tmp1(t,1);
       else
           epsA_skew(t,1) = tmp2(t,1);
       end
    end
    epsA_skew = (epsA_skew - repmat(mean(epsA_skew,1),options_.periods,1))./repmat(std(epsA_skew,1),options_.periods,1); %standardize
    Moments_epsA_skew = zeros(4,1);
    Moments_epsA_skew(1,1) = 0; % as we standardized
    Moments_epsA_skew(2,1) = 1; % as we standardized
    Moments_epsA_skew(3,1) = skewness(epsA_skew); % compute empirical skewness coefficient
    Moments_epsA_skew(4,1) = kurtosis(epsA_skew); % compute empirical kurtosis coefficient
    % use same shocks as in BENCHMARK Gaussian case, but switch out epsA
    idepsA = ismember(M_.exo_names,'epsA');
    exo(:,idepsA) = sqrt(SIGMA2(idepsA,idepsA))*epsA_skew; % scale by the standard deviation    
    SIGMA3(idepsA,idepsA,idepsA) = sqrt(SIGMA2(idepsA,idepsA))^3*Moments_epsA_skew(3,1); % use the third-order product moment and scale it by the standard deviation raised to the 3rd power to get SIGMA3
    oo_ = perturbation_solver_nonsymmetric_order3(M_,oo_,SIGMA2(:),SIGMA3(:));
    RESULTS.(model_variant).exo = exo;
    fprintf('Simulating %s model\n',model_variant)
    RESULTS.(model_variant).y = simult_nonsymmetric_order3(M_,oo_,options_,exo);
end


/*
 **************************************************************************
 * replicate tables
 **************************************************************************
 */

for model_variant = ["BENCHMARK" "CASEI" "CASEII"]
    % Add market price and quantity of risk
    mpr = RESULTS.(model_variant).y(ismember(options_.varobs,'varsdf'),:).*(1+RESULTS.(model_variant).y(ismember(options_.varobs,'R1'),:)/400);
    qrisk = RESULTS.(model_variant).y(ismember(options_.varobs,'TP'),:)./(mpr*400);
    TBL = zeros(length(options_.varobs),4);
    for j = 1:length(options_.varobs)
        TBL(j,1) = mean(RESULTS.(model_variant).y(j,:)');
        TBL(j,2) = std(RESULTS.(model_variant).y(j,:)');
        TBL(j,3) = skewness(RESULTS.(model_variant).y(j,:)');
        TBL(j,4) = kurtosis(RESULTS.(model_variant).y(j,:)');
    end
    TBL = [TBL; mean(RESULTS.(model_variant).exo(:,1)) std(RESULTS.(model_variant).exo(:,1)) skewness(RESULTS.(model_variant).exo(:,1)) kurtosis(RESULTS.(model_variant).exo(:,1))];
    TBL = [TBL; mean(mpr) std(mpr) skewness(mpr) kurtosis(mpr)];
    TBL = [TBL; mean(qrisk) std(qrisk) skewness(qrisk) kurtosis(qrisk)];
    fprintf('    %s\n',upper(model_variant));
    disp(array2table(TBL,'RowNames',[options_.varobs;M_.exo_names(1);'MPR';'Qrisk'],'VariableNames',{'MEAN','STD','SKEWNESS','KURTOSIS'}));
end