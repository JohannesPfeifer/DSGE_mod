/*
 * This file replicates the ABCD-test for the example of the permanent income model provided in
 * Fernandez-Villaverde, Rubio-Ramirez, Sargent, 
 * and Watson (2007), "ABCs (and Ds) of Understanding VARs", American
 * Economic Review, 97(3), 1021-1026
 *
 * Requires the ABCD_test.m to be located in the same folder.
 *
 * This file was written by Johannes Pfeifer.In case you spot mistakes, email me at jpfeifer@gmx.de
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

var y c y_m_c;
    
varexo w;

parameters R sigma_w;

// set parameter values
sigma_w=1;
R= 1.2;

model;
c=c(-1)+sigma_w*(1-R^(-1))*w;
y_m_c=-c(-1)+sigma_w*R^(-1)*w;
y_m_c=y-c;
end;

steady_state_model;
c=0;
y=0;
y_m_c=0;
end;

//set shock variances
shocks;
    var w=1;
end;

steady;
check;
varobs y_m_c; //set to c or y for condition to be satisfied
stoch_simul(order=1,irf=20);

//*************************************perform ABCD-test on model

[result,eigenvalue_modulo,A,B,C,D]=ABCD_test(M_,options_,oo_)

