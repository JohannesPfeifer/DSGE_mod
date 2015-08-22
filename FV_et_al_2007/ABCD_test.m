function [result,eigenvalue_modulo,A,B,C,D]=ABCD_test(M_,options_,oo_)
% function ABCD_test(M_,options_,oo_)
%   computes ABCD Test statistics of Fernandez-Villaverde, Rubio-Ramirez,Sargent, 
%   and Watson (2007), "ABCs (and Ds) of Understanding VARs", American
%   Economic Review, 97(3), 1021-1026
% INPUTS
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   result:             [scalar] dummy indicating whether Conditions 1 (Poor Man's
%                           Invertibility Condition) is satified (1) or not (0)
%   eigenvalue_modulo   [nstates by 1] vector providing the modulo of the eigenvalues of A-B*D^(-1)*C 
%   A                   [nstates by nstates] state transition matrix
%   B                   [nstates by nshocks] matrix with effect of shocks on states 
%   C                   [nobservables by nstates] observation equation mapping matrix of observables into states
%   D                   [nobservables by nshocks] matrix with effect of shocks on observables
%
% Algorithm
%   Uses the state-space system
%           x(+1)=Ax+Bw(+1)
%           y(+1)=Cx+Dw(+1)
%   and computes the eigenvalues of 
%           A-BC^(-1)D  (Condition 1)
%   If all eigenvalues are smaller than 1, the poor man's invertibility
%   condition is satisfied and the structural shocks can be recovered from
%   the observables
% 
%       
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 2013-5 Johannes Pfeifer
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
% You can receive a copy of the GNU General Public License
% at <http://www.gnu.org/licenses/>.

result=NaN;
eigenvalue_modulo=NaN;
A=[];
B=[];
C=[];
D=[];

if ~isfield(options_,'varobs_id')
    warning('ABCD_test: No observables have been defined using the varobs-command.')
    return;    
end

%get state indices
ipred = M_.nstatic+(1:M_.nspred)';
%get state transition matrices
[A,B] = kalman_transition_matrix(oo_.dr,ipred,1:M_.nspred,M_.exo_nbr);
%get observable position in decision rules
obs_var=oo_.dr.inv_order_var(options_.varobs_id);
%get observation equation matrices
[C,D] = kalman_transition_matrix(oo_.dr,obs_var,1:M_.nspred,M_.exo_nbr);
%compute condition


if M_.exo_nbr~=size(options_.varobs_id)
    warning('ABCD_test: ABCD test only works for square case with as many observables as structural shocks.')
    return;
end
if ~isequal(M_.H,0)
    warning('ABCD_test: ABCD test does not work with measurment error not specified as structural shocks.')
    return;
end

if rcond(D)<1e-20
    warning('ABCD_test: Matrix D must not be singular.')
    return;
end
eigenvalues=eig(A-B/D*C);
eigenvalue_modulo=abs(eigenvalues);
if any(eigenvalue_modulo>=1)
    fprintf('\nABCD: The Poor Man''s Invertibility Condition is not satisfied.\n')
    fprintf('ABCD: The structural shocks can NOT be recovered from a VAR in the observables.\n')
    result=0;
else
    fprintf('\nABCD: The Poor Man''s Invertibility Condition is satisfied.\n')
    fprintf('ABCD: The structural shocks can be recovered from a VAR in the observables.\n')
    result=1;
    if any(eigenvalue_modulo>1-10*(eps(eigenvalue_modulo)))
        warning('ABCD: One of the eigenvalues is within machine precision to a unit root. The results may not be accurate.')
    end
end
