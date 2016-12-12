function [shock_mat]=generate_FGRU_shocks(replications,winsorizing_dummy)
% [shock_mat]=generate_FGRU_shocks(replications,winsorizing_dummy)
% generates the FGRU shock matrix using their generation. For Argentina,
% the seed is known and the same as in their study so that the random
% numbers are exactly the same. That in particular implies that the first
% shock vector is the same across draws. The number of periods is
% hardcoded as is the seed.
% 
% 
% Inputs
%   - replications          [scalar] number of replications
%   - winsorizing_dummy     [scalar] dummy indicating whether shocks other than TFP should
%                                       be winsorized
% 
% Outputs:
%   - shock_mat             [96 by 5 by replications] matrix with shocks
% 
% The file is used to replicate the results of
% Benjamin Born and Johannes Pfeifer (2014): "Risk Matters: A comment", American Economic Review
% 
% Copyright (C) 2013-14 Benjamin Born and Johannes Pfeifer
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
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.

randn('seed',2);     %original calibration
shock_1=randn(1,5); %first shock used in FGRU is always the same
temp=randn(95*5*replications,1); %generate other shocks
shock_mat=permute(reshape(temp,5,95,replications),[2 1 3]); %reshape as FGRU filled in row major order
shock_mat=cat(1,repmat(shock_1,[1,1,replications]),shock_mat); %stack first shock that is the same for all draws with other ones
    
%winsorize all shocks except for TFP
if winsorizing_dummy
    winsorizing=(shock_mat>1);
    winsorizing(:,1,:)=0;
    shock_mat(winsorizing)=1;
    winsorizing=(shock_mat<-1);
    winsorizing(:,1,:)=0;
    shock_mat(winsorizing)=-1;
end