function theta_w_opt=get_Calvo_theta(WPC_slope,epsilon_w,betta,varphi,SGU_indicator)
% function theta_w_opt=get_Calvo_theta(WPC_slope,epsilon_w,betta,varphi,SGU_indicator)
% compute the Calvo parameter for a given slope of the wage Phillips Curve
% Inputs:
%   - WPC_slope     [scalar]    slope of the wage Phillips Curve
%   - epsilon_w     [scalar]    substitution elasticity
%   - betta         [scalar]    discount factor
%   - varphi        [scalar]    inverse Frisch elasticity
%   - SGU_indicator [scalar]    1=SGU, 0=EHL

% Copyright (C) 2018 Johannes Pfeifer and Benjamin Born
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
% For a copy of the GNU General Public License,
% see <http://www.gnu.org/licenses/>.

if user_has_matlab_license('optimization_toolbox')
    options=optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
end
if SGU_indicator
    if ~user_has_matlab_license('optimization_toolbox')
        [theta_w_opt,exitflag]=csolve(@(theta_w)WPC_slope-(1-theta_w)*(1-betta*theta_w)/theta_w,0.8,[],1e-10,1000);
        if exitflag==0
            exitflag=1;
        else
            exitflag=-1;
        end       
    else
        [theta_w_opt,fval,exitflag]=fsolve(@(theta_w)WPC_slope-(1-theta_w)*(1-betta*theta_w)/theta_w,0.8,options);
    end
else
    if ~user_has_matlab_license('optimization_toolbox')
        [theta_w_opt,exitflag]=csolve(@(theta_w)WPC_slope-(1-theta_w)*(1-betta*theta_w)/(theta_w*(1+epsilon_w*varphi)),0.8,[],1e-10,1000);
        if exitflag==0
            exitflag=1;
        else
            exitflag=-1;
        end       
    else
        [theta_w_opt,fval,exitflag]=fsolve(@(theta_w)WPC_slope-(1-theta_w)*(1-betta*theta_w)/(theta_w*(1+epsilon_w*varphi)),0.8,options);
    end
end
if theta_w_opt<0 || theta_w_opt>1 || exitflag<1
    error('Could not solve for theta_w')
end