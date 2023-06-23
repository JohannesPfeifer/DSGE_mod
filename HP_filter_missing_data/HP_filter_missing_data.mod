/*
 * This file implements a Hodrick/Prescott HP-filter employing a diffuse Kalman smoother. 
 * It relies on the HP-filter's formulation as a state space problem. See e.g. Section 2 of
 * Hamilton (2018): Why You Should Never Use the Hodrick-Prescott Filter, The Review of 
 * Economics and Statistics 100(5), 831–843.
 *
 * Notes: 
 * - Due to using the Kalman smoother instead of the typical matrix formula, the HP-filter 
 *      naturally handles missing observations/NaN.
 * - The smoothing parameter lambda is given by the ratio of variances var(c)/var(v).
 * - The dataset used contains Argentinean sovereign EMBI spreads from 1993Q4 to 2022Q5, where 
 *      default episodes have been set to NaN/missing.
 *
 * This implementation was written by Johannes Pfeifer. 
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2023 Johannes Pfeifer
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


 
var y   (long_name='observed data') 
    g   (long_name='trend') 
    ;
varexo v (long_name='cyclical component')
    c (long_name='cyclical component')
    ;

parameters lambda (long_name='smoothing parameter')
    ;

lambda=1600;

model(linear);
    [name='observation equation']
    y=g+sqrt(lambda)*c;
    [name='evolution of trend']
    g-g(-1)=g(-1)-g(-2)+v;
end;

steady_state_model;
y=1;
g=1;
end;

steady;

shocks;
    var v=1;
    var c=1;
end;

varobs y;

calib_smoother(datafile=data_spread,diffuse_filter,kalman_algo=1);

temp_data=load('data_spread');
figure
subplot(2,1,1)
plot(oo_.SmoothedVariables.g)
hold on
plot(temp_data.y)
legend({'Trend','Data'})
subplot(2,1,2)
plot(oo_.SmoothedShocks.c)
title('Cyclical component')
