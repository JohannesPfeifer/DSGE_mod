[![DOI](https://zenodo.org/badge/34334463.svg)](https://zenodo.org/badge/latestdoi/34334463)

# DSGE_mod 

A collection of Dynare models. It aims at demonstrating Dynare best practices 
and providing tractable replication files for important models that can be 
useful for further model development. 

## Replicability Issues

The headers of the respective mod-files also note obvious mistakes and typos 
in the respective papers.  

## Compatibility

These mod-files have been tested against Dynare 4.6. Compatibility with 
earlier versions is not guaranteed. For Dynare 4.5 version, please use the 4.5 branch of this repository

# Contributing your own mod-files

Contributions of replication files to this collection are highly welcomed. 
When doing do, e.g. through pull requests, please clearly line out which 
results of the original paper are replicated so that correctness can be 
verified.

# Contained Mod-files

## Aguiar_Gopinath_2007.mod

Replicates Aguiar, Mark and Gopinath, Gita (2007): "Emerging Market Business Cycles: The Cycle is the Trend", Journal of Political Economy, 115(1), pp. 69-102.

This mod-file shows how to deal with trend growth and how to recover the non-stationary variables from the detrended model variables.

## Andreasen_2012

Replicates Table 2 in Andreasen (2012): "On the effects of rare disasters and uncertainty shocks for risk premia in non-linear DSGE models", Review of Economic Dynamics, 15, pp. 295-316.

### Andreasen_2012_rare_disasters.mod
This mod file shows how to simulate DSGE models solved with third-order perturbation and non-symmetric innovations.
The underlying model is a New-Keynesian model with Epstein-Zin preferences.

## Ascari_Sbordone_2014.mod
Replicates Ascari, Guido and Sbordone, Argia M. (2014): "The Macroeconomics of Trend Inflation", Journal of Economic Literature, 52(3), pp. 679-739.

This mod-file shows how to access steady state variables in order to plot steady state dependencies on parameters. 
It also shows how to manually do a stability mapping be iterating over a grid on the parameter space.
 
## Basu_Bundick_2017.mod

Replicates the Generalized Impulse Response Functions (GIRFs) at the stochastic steady/ergodic mean in the absence of shocks by Basu/Bundick (2017): 
"Uncertainty shocks in a model of effective demand", Econometrica, 85(3), pp. 937-958

## Born_Pfeifer_2014

Replicates Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", 
American Economic Review, 104(12), pp. 4231-4239.

### Born_Pfeifer_RM_Comment.mod

This mod-file shows how to estimate a model solved with third order perturbation using the
Simulated Method of Moments. It also shows how to generate IRFs at the stochastic steady state (ergodic
mean in the absence of shocks (EMAS) in the terminology of the paper). For practical purposes it is
highly recommended to use the standard Andreasen et al. (2013) pruning scheme available in Dynare's `simult_.m` 
instead of the FGRU version in `simult_FGRU.m` (see the comments in the mod-file).

## Born_Pfeifer_2018

Replicates Benjamin Born and Johannes Pfeifer (2018): "The New Keynesian Wage Phillips Curve: Calvo vs. Rotemberg",
Macroeconomic Dynamics, 24, 2020, 1017–1041.

### Monetary_Policy_IRFs/Born_Pfeifer_2018_MP
`run_IRF_comparison.m` creates "Figure 1: Impulse response functions to 1 percentage point 
(annualized) monetary policy shock under Calvo".

### Welfare/Born_Pfeifer_2018_welfare.mod
The mod-file `Born_Pfeifer_2018_welfare.mod` shows how to compute conditional and unconditional welfare.
`run_welfare_comparison_efficient_steady_state.m` and `run_welfare_comparison_inefficient_steady_state.m` create
the welfare comparison between the four different labor market setups presented in Tables 4 and 5 of the paper.

## Born_Pfeifer_2020

Replicates the DSGE model results of Benjamin Born and Johannes Pfeifer (2020): "Uncertainty-driven business cycles:
assessing the markup channel", forthcoming at Quantitative Economics. The main file is run_model_IRF_generation.m

Note that sequential calling of Dynare can cause problems on Windows if the created files are temporarily locked by other 
processes like e.g. cloud drive apps. We recommend not running the codes in folders synchronized by cloud drives.

## Caldara_et_al_2012.mod

Replicates Caldara, Dario and Fernandez-Villaverde, Jesus and Rubio-Ramirez, 
Juan F. and Yao, Wen (2012): "Computing DSGE Models with Recursive 
Preferences and Stochastic Volatility", Review of Economic Dynamics, 15, pp. 
188-206. 

This mod-file shows how to use auxiliary variables to deal with recursive preferences and expected returns. 
It also shows how to use the ```plot_policy_fun.m``` to plot the policy functions using Dynare

## Chari_et_al_2007.mod

Replicates Chari, V.V/Kehoe, Patrick J./McGrattan, Ellen (2007), "Business Cycle Accounting", Econometrica, 75(3), pp. 781-836. 

It demonstrates how to use the linearized benchmark model estimated using Maximum Likelihood to conduct the 
Business Cycle Accounting  as is done in the paper for the 1982 recession.

## Collard_2001

### Collard_2001_example1.mod
Replicates example 1 of Collard (2001): "Stochastic simulations with DYNARE: A practical guide".
The file `get_shock_standard_deviation` shows how to find the shock size that induces a given response 
of an endogenous variable in a specified period.

## FV_et_al_2007 

Provides codes for the ABCD-test of Fernandez-Villaverde, Rubio-Ramirez, Sargent, and Watson (2007), 
"ABCs (and Ds) of Understanding VARs", American Economic Review, 97(3), pp. 
1021-1026 

Includes the ```ABCD_test.m```. Note that it tests only a sufficient condition, not
a necessary one, if the minimal state space is not computed. For details, see e.g. Komunjer/Ng (2011):
"Dynamic Identification of Dynamic Stochastic General Equilibrium Models", Econometrica, 79(6), 1995–2032.

### FV_et_al_2007.mod

Replicates the ABCD-test for the example of the permanent income model 
provided in Fernandez-Villaverde, Rubio-Ramirez, Sargent, and Watson (2007), 
"ABCs (and Ds) of Understanding VARs", American Economic Review, 97(3), pp. 
1021-1026

### FV_et_al_2007_ABCD_minreal.mod

Shows how to compute the minimal state space using Matlab's Control toolbox for the example of 
Saccal, Alessandro (2020): "A note on minimality in Dynare", available at https://mpra.ub.uni-muenchen.de/103656/1/MPRA_paper_103656.pdf.

## Gali_2008

### Gali_2008_chapter_2.mod

Implements the baseline Classical Monetary Economy model of Jordi Galí­ 
(2008): Monetary Policy, Inflation, and the Business Cycle, Princeton 
University Press, Chapter 2 

### Gali_2008_chapter_3.mod 

Implements the baseline New Keynesian model of Jordi Galí (2008): Monetary  
Policy, Inflation, and the Business Cycle, Princeton University Press, 
Chapter 3 

### Gali_2008_chapter_5_discretion.mod 

Implements the optimal monetary policy under discretion exercise of Jordi 
Galí (2008): Monetary  Policy, Inflation, and the Business Cycle, Princeton 
University Press, Chapter 5.1.1. It shows how to use the 
```discretionary_policy``` command.

### Gali_2008_chapter_5_commitment.mod 

Implements the optimal monetary policy under commitment exercise of Jordi 
Galí (2008): Monetary  Policy, Inflation, and the Business Cycle, Princeton 
University Press, Chapter 5.1.2. It shows how to use the ```ramsey_policy``` 
command. 

## Gali_2010

### Gali_2010.mod 
This file was written together with Lahcen Bounader. It replicates the results of the baseline sticky wage model 
of Jordi Galí (2010): Monetary Policy and Unemployment, Handbook of Monetary Economics, Volume 3A, 
Chapter 10, pp. 487-546.
Please see the header of the mod-file for additional remarks.

### Gali_2010_calib_target.mod 
This file was written together with Lahcen Bounader. It implements the baseline sticky wage model 
of Jordi Galí (2010): Monetary Policy and Unemployment, Handbook of Monetary Economics, Volume 3A, 
Chapter 10, pp. 487-546. When doing so, it corrects issues with the original calibration of Gali (2010).
It demonstrates how in a linearized model a steady_state-file can be used to set the deep parameters of the
model to satisfy calibration targets on the non-linear model. The steady_state-file takes the calibration targets 
and calls a numerical solver on some of the nonlinear steady state equations to get the corresponding parameters 
that make the steady state satisfy the targets.

## Gali_2015

### Gali_2015_chapter_2.mod 
Implements the baseline Classical Monetary Economy of Jordi Gali (2015): Monetary Policy, 
Inflation, and the Business Cycle, Princeton University Press, Second Edition, Chapter 2

### Gali_2015_chapter_3.mod and Gali_2015_chapter_3_nonlinear.mod 
Implements the baseline New Keynesian model of Jordi Gali (2015): Monetary Policy, 
Inflation, and the Business Cycle, Princeton University Press, Second Edition, Chapter 3 

### Gali_2015_chapter_4.mod
Implements the welfare analysis of Chapter 4.4 on simple rules in the baseline New Keynesian model
of Jordi Gali (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press, Second Edition.

### Gali_2015_chapter_5_commitment.mod 

Implements the optimal monetary policy under commitment exercise of Jordi Gali (2015): Monetary Policy, Inflation, 
and the Business Cycle, Princeton University Press, Second Edition, Chapter 5.2.2. It shows how to use the ```ramsey_policy``` 
command. 

### Gali_2015_chapter_5_commitment_ZLB.mod 

Implements the optimal monetary policy at the ZLB under commitment exercise 
Jordi Gali (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press, Second Edition,  
Chapter 5.4.2. It shows how to solve a perfect foresight model with a Levenberg-Marquardt mixed complementarity problem (lmmcp)
approach to deal with the zero lower bound on interest rates.

### Gali_2015_chapter_5_discretion.mod 

Implements the optimal monetary policy under discretion exercise of Jordi Gali (2015): Monetary Policy, Inflation, 
and the Business Cycle, Princeton University Press, Second Edition, Chapter 5.2.1. It shows how to use the 
```discretionary_policy``` command.

### Gali_2015_chapter_5_discretion_ZLB.mod 

Implements the optimal monetary policy at the ZLB under discretion exercise 
Jordi Gali (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press, Second Edition,  
Chapter 5.4.1. It shows how to solve a perfect foresight model with a Levenberg-Marquardt mixed complementarity problem (lmmcp)
approach to deal with the zero lower bound on interest rates.

### Gali_2015_chapter_6.mod
Implements the New Keynesian model with price and wage rigidities 
of Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton 
University Press, Second Edition, Chapter 6

### Gali_2015_chapter_6_4.mod
Implements the New Keynesian model with price and wage rigidities under optimal policy 
with commitment (Ramsey) of Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton 
University Press, Second Edition, Chapter 6.4

### Gali_2015_chapter_6_5.mod
Implements the New Keynesian model with price and wage rigidities under under simple rules
of Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton 
University Press, Second Edition, Chapter 6.5

### Gali_2015_chapter_7.mod
Implements the New Keynesian model with price and wage rigidities and unemployment 
of Chapter 7 of Jordi Galí (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton 
University Press, Second Edition.

### Gali_2015_chapter_8.mod
Implements the baseline New Keynesian small open economy model of Chapter 8 of
Jordi Gali (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press, Second Edition.

## Gali_Monacelli_2005.mod 

Replicates Galí, Jordi and Monacelli, Tommaso (2005): "Monetary Policy and 
Exchange Rate Volatility in a Small Open Economy", Review of Economic Studies 
72, pp. 707-734. 

This mod-file shows how to use Dynare's LaTeX-capacities

## GarciaCicco_et_al_2010.mod

Replicates the model studied in García-Cicco, Javier and Pancrazi, Roberto and Uribe, Martín (2010): 
"Real Business Cycles in Emerging Countries", American Economic Review, 100(5), pp. 2510-2531.
 
It provides a replication code for the main results of the original paper for the case of Argentina. 

This mod-file shows how to use the loglinear and logdata options of Dynare.

## Ghironi_Melitz_2005.mod

This file replicates the Baseline model under Financial Autarky of
Ghironi/Melitz (2005), "International trade and macroeconomic dynamics with 
heterogeneous firms", Quarterly Journal of Economics, 120(3), 865-915.

## Guerrieri_Iacoviello_2015
Replicates the example results of Guerrieri, Luca and Iacoviello, Matteo (2015): 
"OccBin: A toolkit for solving dynamic models with occastionally binding
constraints easily", Journal of Monetary Economics 70, pp.22-38.
It provides replication codes for the IRFs in figures 3 and 5.
The mod files show how to use Dynare's occbin toolbox for stochastic 
simulations with occasionally binding constraints.

### Guerrieri_Iacoviello_2015_rbc.mod
Replicates the RBC model with a constraint on investment (irreversible investment).

### Guerrieri_Iacoviello_2015_nk.mod
Replicates the New Keynesian model with a constraint on the nominal interest rate (zero-lower-bound).

## HP_filter_missing_data.mod

This file implements a Hodrick/Prescott HP-filter employing a diffuse Kalman smoother. Due 
to using the Kalman smoother instead of the typical matrix formula, the HP-filter 
naturally handles missing observations/NaN.

## Hansen_1985.mod

Replicates the model studied in Hansen, Gary D. (1985): "Invisible labor and 
the business cycle", Journal of Monetary Economics 16, pp.309-327.

This mod-file shows how to use the loglinear option to get moments of 
percentage deviations without loglinearizing the model and how to use the ```
get_simul_replications.m``` file to read out simulations generated by the 
```simul_replic``` option 

## Ireland_2004.mod

Estimates the New Keynesian model of Ireland, Peter (2004): "Technology shocks in the New Keynesian Model",
Review of Economics and Statistics, 86(4), pp. 923-936

This mod-file shows how to estimate DSGE models using maximum likelihood in Dynare.

## Jermann_1998.mod

This file replicates some of the results in Jermann (1998): Asset pricing in production economies,
Journal of Monetary Economics, 41, pp. 257-275, using a second-order perturbation approximation.

## Jermann_Quadrini_2012

Provides replication files for Pfeifer, Johannes (2016): "Macroeconomic effects effects of financial 
shocks: A comment", Dynare Working Paper 50. This paper replicates and corrects the results obtained in 
Jermann/Quadrini (2012): "Macroeconomic effects of financial shocks", American Economic Review, 102(1): 
238-271.

### Jermann_Quadrini_2012_RBC
Implements the RBC model of Jermann/Quadrini (2012). It allows replicating the original results and 
generates the results of Pfeifer (2016), who documents a mistake in the TFP-construction of JQ 
that requires recalibrating the model.

### Jermann_Quadrini_2012_NK
This file replicates the estimation of the New Keynesian model of Jermann/Quadrini (2012) conducted and 
described in Pfeifer (2016).

## McCandless_2008

This folder contains replication files for George McCandless (2008): The ABCs of RBCs - An Introduction to Dynamic 
Macroeconomic Models, Harvard University Press

### McCandless_2008_Chapter_9.mod

This file replicates the Money in the Utility Function model studied in Chapter 9

### McCandless_2008_Chapter_9.mod

This file replicates the open economy model studied in Chapter 13

## NK_linear_forward_guidance.mod

Shows how to implement forward guidance in a baseline New Keynesian model using a sequence of monetary policy shocks.

## RBC_IRF_matching.mod

This file takes the baseline RBC model with TFP and government spending shocks, calibrated to US data from 1947Q4:2016Q1 and
estimates the persistence of the AR(2) government spending shock via impulse response function (IRF) matching. 

## RBC_baseline

### RBC_baseline.mod

This file presents a baseline RBC model with TFP and government spending shocks, calibrated to US data from 1947Q4:2016Q1. 
The model setup is described in Handout_RBC_model.pdf and resembles the one in King/Rebelo (1999): Resuscitating Real Business Cycles, Handbook of Macroeconomics, Volume 1, and
Romer (2012), Advanced macroeconomics, 4th edition. The driving processes are estimated as AR(1)-processes on linearly detrended data.

### RBC_baseline_first_diff_bayesian.mod

Estimates the baseline RBC model on simulated data.

## RBC_baseline_welfare.mod

Computes the welfare-maximizing optimal labor tax rate in a baseline RBC model with only TFP shocks. It does
so by defining welfare recursively in the model block and calling an optimizer to find the parameter for the 
steady state tax rate that maximizes welfare.

Estimates the baseline RBC model on simulated data.

## RBC_capitalstock_shock.mod

Implements a simple RBC model with a time t shock to the capital stock. 

## RBC_news_shock_model.mod

Implements a simple RBC model with additively separable utility and TFP news 
calibrated to US data. It shows how to generate IRFs to a "pure" news shock 
where an 8 period anticipated news shock does not materialize at time 0. This 
is the type of policy experiment that is for example performed in Beaudry 
Portier (2004): An exploration into Pigou's theory of cycles, Journal of 
Monetary Economics 51, pp. 1183-1216. 

## RBC_state_dependent_GIRF.mod

This file takes the baseline RBC model and demonstrates how to compute Generalized 
Impulse Response Functions using Dynare's simult_-function. The model is solved up 
to second order to allow for non-linearities.

## Ramsey_Cass_Koopmans.mod 
Studies the transition behavior of a simple Ramsey/Cass/Koopmans economy with Cobb-Douglas production 
function to its balanced growth path (BGP). The RCK model is solved here in aggregate, 
i.e. non-detrended form along its balanced growth path. For that purpose, trending labor-augmenting
technology and population processes are defined.

## SGU_2003.mod 

Replicates Schmitt-Grohé, Stephanie and Uribe, Martín (2003): "Closing small 
open economy models", Journal of International Economics, 61, pp. 163-185. 

## SGU_2004.mod

Replicates the neoclassical growth model for Schmitt-Grohé/Uribe (2004): 
"Solving dynamic general equilibrium models using a second-order 
approximation to the policy function", Journal of Economic Dynamics & 
Control, 28, pp. 755-775 

## Sims_2012

### Sims_2012_RBC.mod

Replicates the results for the basic RBC model presented in Eric R. Sims 
(2012): "New, Non-Invertibility, and Structural VARs", Advances in 
Econometrics, Volume 28, 81-135 

Requires the ABCD_test.m from the FV_et_al_2007-folder to be located in the 
same folder. 

## Smets Wouters 2007

Provides replication files for Smets, Frank and Wouters, Rafael (2007):  
"Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach", 
American Economic Review, 97(3), pp. 586-606.

### Smets_Wouters_2007.mod 

Rudimentary code that is compatible with Dynare 4.2.5 onwards. See also the
header to Smets_Wouters_2007_45.mod for additional remarks.

### Smets_Wouters_2007_45.mod 

Provides replication files that are compatible with Dynare 
4.5 onwards and make full use of Dynare's LaTeX-capabilities to better 
document the original replication files. 

## Solow_model

Various mod-files related to the basic Solow-Swan model, using Dynare's perfect
foresight routines to study steady state transitions when e.g. parameters change

### Solow_SS_transition.mod 
Studies the transition behavior of a simple Solow-Swan economy with Cobb-Douglas 
production function to its steady state when started with a capital stock different
from steady state.

### Solow_growth_rate_changes.mod 
Studies the transition behavior of a simple Solow-Swan economy with Cobb-Douglas production 
function after unanticipated for changes in technology or population growth. 

### Solow_nonstationary.mod 
Studies the transition behavior of a simple Solow-Swan economy with Cobb-Douglas production 
function to its balanced growth path (BGP). The Solow model is solved here in aggregate, 
i.e. non-detrended form along its balanced growth path. For that purpose, trending labor-augmenting
technology and population processes are defined.

## Stock_SIR_2020.mod

This file implements a simple Susceptible-Infected-Recovered (SIR) model as in 
James H. Stock (2020): "Data Gaps and the Policy Response to the Novel Coronavirus". 

## Woodford_2003

### Woodford_2003_Chapter_7.mod

Implements the deterministic optimal policy exercise in Figure 7.1 of Michael Woodford (2003): 
"Interest and prices", Princeton University Press, page 476. The same figure is reproduced as 
Figure 2 in Michael Woodford (2010): "Optimal Monetary Stabilization Policy", Chapter 14 of the 
Handbook of Monetary Economics, Volume 3B, Elsevier
