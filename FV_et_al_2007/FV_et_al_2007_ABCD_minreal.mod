/*
 * This file shows that the original ABCD-test is too strict as it does not check minimality (see 
 * e.g. Komunjer/Ng (2011): "Dynamic Identification of stochastic general equilibrium models", 
 * Econometrica, 79(6), 1995--2032). It employs a different version of the permanent income model of 
 * Fernandez-Villaverde, Rubio-Ramirez, Sargent, and Watson (2007), "ABCs (and Ds) of Understanding VARs", 
 * American Economic Review, 97(3), 1021-1026, based on Saccal, Alessandro (2020): "A note on minimality in Dynare", 
 * available at https://mpra.ub.uni-muenchen.de/103656/1/MPRA_paper_103656.pdf.
 * In this example, the original ABCD-test will fail, because the relevant state-space system is not minimal: 
 * The observed y_pt has a 0 coefficient on c_pt(-1) (in contrast to c, which is why it is a state for the 
 * full system). Thus, the state is immaterial. 
 * 
 * Using the minimal system based on Matlab's minreal instead returns the 
 * correct result. The example is  
 * 
 *
 * Requires the ABCD_test.m to be located in the same folder.
 *
 * This file was written by Johannes Pfeifer. In case you spot mistakes, email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright (C) 2013 Johannes Pfeifer
 *
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can receive a copy of the GNU General Public License
 * at <http://www.gnu.org/licenses/>.
 */

var c y_pt s;
    
varexo w;

parameters R sigma_w;

// set parameter values
sigma_w=1;
R= 1.2;

model;
c=c(-1)+sigma_w*(1-R^(-1))*w;
y_pt=sigma_w*w;
s=y_pt-c;
end;

steady_state_model;
c=0;
y_pt=0;
s=0;
end;

//set shock variances
shocks;
    var w=1;
end;

steady;
check;
varobs y_pt; //set to c or y for condition to be satisfied
stoch_simul(order=1,irf=20);

//*************************************perform ABCD-test on model

[result,eigenvalue_modulo,A,B,C,D]=ABCD_test(M_,options_,oo_,0); 
[result,eigenvalue_modulo,A,B,C,D]=ABCD_test(M_,options_,oo_,1); %check minimality

