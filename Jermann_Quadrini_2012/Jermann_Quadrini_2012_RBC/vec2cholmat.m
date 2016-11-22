function L_chol=vec2cholmat(chol_vec)
%L_chol=vec2cholmat(chol_vec)
%Writes a vector containing the lower diagonal elements of a Cholesky
%decomposed matrix back into a matrix
%
% Inputs:
%   - chol_vec [nvars*(nvars-1)/2] vector containing the lower diagonal of
%           a Cholesky decomposed matrix; elements are in column order, i.e. (1,1),
%           (2,1), ..., (1,2), (2,2),...
% Outputs:
%   - L_chol    [nvars by nvars] Lower triangular Cholesky matrix

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

nvars=round(-0.5+sqrt(1/4+length(chol_vec)*2));

L_chol=zeros(nvars,nvars);
iter=1;
for column=1:nvars
    for row=column:nvars
        L_chol(row,column)=chol_vec(iter,1);
        iter=iter+1;
    end
end
