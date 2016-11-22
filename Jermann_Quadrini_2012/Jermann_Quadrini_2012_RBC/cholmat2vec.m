function chol_vec=cholmat2vec(Q_chol)
%chol_vec=vec2cholmat(chol_vec)
%Writes the lower diagonal elements of a Cholesky
%decomposed matrix back into a vector
%
% Inputs:
%   - Q_chol    [nvars by nvars] Lower triangular Cholesky matrix
%
% Outputs:
%   - chol_vec [nvars*(nvars-1)/2] vector containing the lower diagonal of
%           a Cholesky decomposed matrix; elements are in column order, i.e. (1,1),
%           (2,1), ..., (1,2), (2,2),...
% 
% Copyright (C) 2014-2016 Johannes Pfeifer
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

nstates=size(Q_chol,1);
chol_vec=zeros(nstates*(nstates+1)/2,1);
iter=1;
for column=1:nstates
    for row=column:nstates
        chol_vec(iter,1)=Q_chol(row,column);
        iter=iter+1;
    end
end
