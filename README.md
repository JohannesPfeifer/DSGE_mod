# DSGE_mod 

A collection of Dynare models. It aims at demonstrating Dynare best practices 
and providing tractable replication files for important models that can be 
useful for further model development. 

## Replicability Issues

The headers of the respective mod-files also note obvious mistakes and typos 
in the respective papers.  

## Compatibility

These mod-files have been tested against the current unstable version, to be released
as Dynare 4.5. Compatibility with earlier versions is not guaranteed.

# Contributing your own mod-files

Contributions of replication files to this collection are highly welcomed. 
When doing do, e.g. through pull requests, please clearly line out which 
results of the original paper are replicated so that correctness can be 
verified.

# Contained Mod-files

## Ascari_Sbordone_2014.mod
Replicates Ascari, Guido and Sbordone, Argia M. (2014): "The Macroeconomics of Trend Inflation", Journal of Economic Literature, 52(3), pp. 679-739.

This mod-file shows how to access steady state variables in order to plot steady state dependences on parameters. 
It also shows how to manually do a stability mapping be iterating over a grid on the parameter space.
 
## Aguiar_Gopinath_2007.mod

Replicates Aguiar, Mark and Gopinath, Gita (2007): "Emerging Market Business Cycles: The Cycle is the Trend", Journal of Political Economy, 115(1), pp. 69-102.

This mod-file shows how to deal with trend growth and how to recover the non-stationary variables from the detrended model variables.

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

## FV_et_al_2007 

Replicates the ABCD-test for the example of the permanent income model 
provided in Fernandez-Villaverde, Rubio-Ramirez, Sargent, and Watson (2007), 
"ABCs (and Ds) of Understanding VARs", American Economic Review, 97(3), pp. 
1021-1026 

Includes the ```ABCD_test.m```

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

## Gali_2015

### Gali_2015_chapter_4.mod
Implements the welfare analysis of Chapter 4.4 on simple rules in the baseline New Keynesian model
of Jordi Gali (2015): Monetary Policy, Inflation, and the Business Cycle, Princeton University Press, Second Edition.

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

## RBC_capitalstock_shock.mod

Implements a simple RBC model with a time t shock to the capital stock. 


## RBC_news_shock_model.mod

Implements a simple RBC model with additively separable utility and TFP news 
calibrated to US data. It shows how to generate IRFs to a "pure" news shock 
where an 8 period anticipated news shock does not materialize at time 0. This 
is the type of policy experiment that is for example performed in Beaudry 
Portier (2004): An exploration into Pigou's theory of cycles, Journal of 
Monetary Economics 51, pp. 1183-1216. 

## SGU_2003.mod 

Replicates Schmitt-Grohé, Stephanie and Uribe, Martín (2003): "Closing small 
open economy models", Journal of International Economics, 61, pp. 163-185. 

## SGU_2004.mod

Replicates the neoclassical growth model for Schmitt-Grohé/Uribe (2004): 
"Solving dynamic general equilibrium models using a second-order 
approximation to the policy function", Journal of Economic Dynamics & 
Control, 28, pp. 755-775 

## Sims_2012 

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

## Smets_Wouters_2007_45.mod 

Provides replication files that are compatible with Dynare 
4.5 onwards and make full use of Dynare's LaTeX-capabilities to better 
document the original replication files. 
