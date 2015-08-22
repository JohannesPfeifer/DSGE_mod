# DSGE_mod 

A collection of Dynare models. It aims at demonstrating Dynare best practices 
and providing tractable replication files for important models that can be 
useful for further model development. 

# Contained Mod-files

## Ascari_Sbordone_2014.mod
Replicates Ascari, Guido and Sbordone, Argia M. (2014): "The Macroeconomics of Trend Inflation", Journal of Economic Literature, 52(3), pp. 679-739.

This mod-file shows how to access steady state variables in order to plot steady state dependences on parameters. It also shows how to manually do a stability mapping be iterating over a grid on the parameter space.
 
## Aguiar_Gopinath_2007.mod

Replicates Aguiar, Mark and Gopinath, Gita (2007): "Emerging Market Business Cycles: The Cycle is the Trend", Journal of Political Economy, 115(1), pp. 69-102.

This mod-file shows how to deal with trend growth and how to recover the non-stationary variables from the detrended model variables.

## Caldara_et_al_2012.mod

Replicates Caldara, Dario and Fernandez-Villaverde, Jesus and Rubio-Ramirez, 
Juan F. and Yao, Wen (2012): "Computing DSGE Models with Recursive 
Preferences and Stochastic Volatility", Review of Economic Dynamics, 15, pp. 
188-206. 

This mod-file shows how to use auxiliary variables to deal with recursive preferences and expected returns. It also shows how to use the ```plot_policy_fun.m``` to plot the policy functions using Dynare

## FV_et_al_2007 

Replicates the ABCD-test for the example of the permanent income model 
provided in Fernandez-Villaverde, Rubio-Ramirez, Sargent, and Watson (2007), 
"ABCs (and Ds) of Understanding VARs", American Economic Review, 97(3), pp. 
1021-1026 

Includes the ```ABCD_test.m```

## Gali_2008

### Gali_2008_chapter_2.mod

Implements the baseline Classical Monetary Economy model of Jordi Galí 
(2008): Monetary Policy, Inflation, and the Business Cycle, Princeton 
University Press, Chapter 2 

### Gali_2008_chapter_3.mod 

Implements the baseline New Keynesian model of Jordi Galí (2008): Monetary  
Policy, Inflation, and the Business Cycle, Princeton University Press, 
Chapter 3 

## Gali_Monacelli_2014.mod 

Replicates Galí, Jordi and Monacelli, Tommaso (2005): "Monetary Policy and 
Exchange Rate Volatility in a Small Open Economy", Review of Economic Studies 
72, pp. 707-734. 

This mod-file shows how to use Dynare's LaTeX-capacities

## RBC_capitalstock_shock.mod

Implements a simple RBC model with a time t shock to the capital stock. 


## RBC_news_shock_model.mod

Implements a simple RBC model with additively separable utility and TFP news 
calibrated to US data. It shows how to generate IRFs to a "pure" news shock 
where an 8 period anticipated news shock does not materialize at time 0. This 
is the type of policy experiment that is for example performed in Beaudry 
Portier (2004): An exploration into Pigou’s theory of cycles, Journal of 
Monetary Economics 51, pp. 1183–1216. 

## SGU_2003.mod 

Replicates Schmitt-Grohé, Stephanie and Uribe, Martín (2003): "Closing small 
open economy models", Journal of International Economics, 61, pp. 163-185. 

## SGU_2004.mod

Replicates the neoclassical growth model for Schmitt-Grohé/Uribe (2004): 
"Solving dynamic general equilibrium models using a second-order 
approximation to the policy function", Journal of Economic Dynamics & 
Control, 28, pp. 755 – 775 

## Sims_2012 

Replicates the results for the basic RBC model presented in Eric R. Sims 
(2012): "New, Non-Invertibility, and Structural VARs", Advances in 
Econometrics, Volume 28, 81–135 

Requires the ABCD_test.m from the FV_et_al_2007-folder to be located in the 
same folder. 

## Smets_Wouters_2007.mod 

Provides replication files for Smets, Frank and Wouters, Rafael (2007):  
"Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach", 
American Economic Review, 97(3), pp. 586-606, that are compatible with Dynare 
4.2.5 onwards 
