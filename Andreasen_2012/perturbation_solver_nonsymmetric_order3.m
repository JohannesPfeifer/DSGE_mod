function oo_ = perturbation_solver_nonsymmetric_order3(M_,oo_,SIGMA2,SIGMA3)
% oo_ = perturbation_solver_nonsymmetric_order3(M_,oo_,SIGMA2,SIGMA3)
% -------------------------------------------------------------------------
% This function recomputes the perturbation matrices ghs2, ghs3, ghxss, ghuss
% in case of non-Gaussian innovations specified by nonzero third-order 
% product moments SIGMA3 = E[kron(u,kron(u,u))] and possibly second-order
% product moments SIGMA2 = E[kron(u,u)] that are different to what is declared
% in the shocks block.
% -------------------------------------------------------------------------
% NOTES
% - Dynare assumes Gaussian innovations; hence, SIGMA3=0 and ghs3=0 by default.
%   Therefore if you want to do stochastic simulations or irfs, you have to
%   write your own codes and not use Dynare's simult_.m function, see for
%   example simult_nonsymmetric_order3.m
% - ghs2, ghxss, ghuss need to be re-computed only if SIGMA2≠M_.Sigma_e(:)
% - this implementation is not really efficient; one should probably exploitmake
%   sparsity and symmetry by e.g. using sparse_hessian_times_B_kronecker_C and A_times_B_kronecker_C
% NOTATION
% - y0 denotes CURRENT endogenous variables at t
% - y_ denotes PREVIOUS state variables (predetermined and mixed) at t-1
% - yp denotes FUTURE jumper variables (mixed and forward) at t+1
% - u0 denotes current shocks
% - up denotes future shocks
% - s denotes perturbation parameter
% - z is vector of dynamic variables + exogenous variables, i.e. z=[y_' y0' yp' u0']',
%   used to evaluate the partial derivatives of the dynamic model
% =========================================================================
% INPUTS
%  * M_:                     [structure] information about model
%  * oo_:                    [structure] storage for results (particularly lower-order perturbation matrices)
%  * SIGMA2:                 [M_.exo_nbr^2 x 1] vector of second-order product moments of shocks
%  * SIGMA3:                 [M_.exo_nbr^3 x 1] vector of third-order product moments of shocks
% -------------------------------------------------------------------------
% OUTPUTS
%  * oo_:                    [structure] storage for updated results
% -------------------------------------------------------------------------
% THIS FUNCTION CALLS
%   * [fname,'.dynamic']
%   * gensylv
%   * unfold_g3
% =========================================================================
% Copyright (C) 2022 Willi Mutschler (@wmutschl, willi@mutschler.eu)
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
% see <https://www.gnu.org/licenses/>.
% =========================================================================

%% some useful indices
[I,~] = find(M_.lead_lag_incidence');              % index for dynamic variables that actually appear in the model, in declaration order
y_ = transpose(M_.nstatic+(1:M_.nspred));          % index for state variables in DR order
yp = transpose(M_.nstatic+M_.npred+(1:M_.nsfwrd)); % index for jumper variables in DR order
z_nbr = nnz(M_.lead_lag_incidence) + M_.exo_nbr;   % number of dynamic variables in Jacobian (including exogenous)
z1 = [nonzeros(M_.lead_lag_incidence(:,oo_.dr.order_var)'); nnz(M_.lead_lag_incidence)+(1:M_.exo_nbr)']; % index for columns in first dynamic derivative matrix
z2 = permute(reshape(1:z_nbr^2,z_nbr,z_nbr),[2 1]); % index for columns in second dynamic derivative matrix g2
z3 = permute(reshape(1:z_nbr^3,z_nbr,z_nbr,z_nbr),[3 2 1]); % index for columns in third dynamic derivative matrix g3

%% dynamic Jacobians
[~, d1f, d2f, d3f] = feval([M_.fname,'.dynamic'], oo_.steady_state(I), oo_.exo_steady_state', M_.params, oo_.steady_state, 1); % get derivatives of dynamic model from script files evaluated at non-stochastic steady-state in declaration order
d3f = identification.unfold_g3(d3f,z_nbr); % d3f does not contain symmetric elements, so we need to unfold it
d1f_d1y0 = d1f(:,nonzeros(M_.lead_lag_incidence(2,oo_.dr.order_var))); % first derivative of dynamic model with respect to y0
d1f_d1yp = d1f(:,nonzeros(M_.lead_lag_incidence(3,oo_.dr.order_var))); % first derivative of dynamic model with respect to yp (one time)
d2f_d2yp = d2f(:,z2(nonzeros(M_.lead_lag_incidence(3,oo_.dr.order_var)),nonzeros(M_.lead_lag_incidence(3,oo_.dr.order_var)))); % second derivative of dynamic model with respect to yp (two times)
d3f_d3yp = d3f(:,z3(nonzeros(M_.lead_lag_incidence(3,oo_.dr.order_var)),nonzeros(M_.lead_lag_incidence(3,oo_.dr.order_var)),nonzeros(M_.lead_lag_incidence(3,oo_.dr.order_var)))); % third derivative of dynamic model with respect to yp (three times)
d2f = d2f(:,z2(z1,z1)); % second partial derivative matrix in DR order
d3f = d3f(:,z3(z1,z1,z1)); % third partial derivative matrix in DR order

%% canoncial perturbation matrices that are needed at higher orders
A = zeros(M_.endo_nbr,M_.endo_nbr);
A(:,M_.lead_lag_incidence(2,oo_.dr.order_var)~=0) = d1f_d1y0;
A(:,M_.nstatic+(1:M_.nspred)) = A(:,M_.nstatic+(1:M_.nspred)) + d1f_d1yp*oo_.dr.ghx(yp,:);
B = zeros(M_.endo_nbr,M_.endo_nbr);
B(:,M_.nstatic+M_.npred+1:end) = d1f_d1yp;

%% Re-compute ghs2
if ~isequal(SIGMA2,M_.Sigma_e(:))
    RHS = (d1f_d1yp*oo_.dr.ghuu(yp,:)+d2f_d2yp*kron(oo_.dr.ghu(yp,:),oo_.dr.ghu(yp,:)))*SIGMA2;
    ghs2 = -(A+B)\RHS;    
end

%% Re-compute ghxss
if ~isequal(SIGMA2,M_.Sigma_e(:))
    % permuted identity matrices useful for tensor unfolding, see e.g. Levintal (2017)
    I_uxu = speye(M_.exo_nbr * M_.nspred  * M_.exo_nbr );
    id_uxu = reshape(1:(M_.exo_nbr * M_.nspred  * M_.exo_nbr ),1,  M_.exo_nbr , M_.nspred  , M_.exo_nbr );
    P_up1_x1up2 = I_uxu(:,ipermute(id_uxu,[0 flip(-[2,1,3]+4)]+1)); %this permutes e.g. kron(zup,zxup) to make it compatible with [zup]_d1[zxup]_a1d2 tensor
    P_up2_x1up1 = I_uxu(:,ipermute(id_uxu,[0 flip(-[3,1,2]+4)]+1)); %this permutes e.g. kron(zup,zxup) to make it compatible with [zup]_d2[zxup]_a1d1 tensor

    Ix = eye(M_.nspred);
    Iu = eye(M_.exo_nbr);
    Iu_kron_Iu = eye(M_.exo_nbr^2);
    zx = [Ix;
          oo_.dr.ghx;
          oo_.dr.ghx(yp,:)*oo_.dr.ghx(y_,:);
          zeros(M_.exo_nbr,M_.nspred)];
    zup = [zeros(M_.nspred,M_.exo_nbr);
           zeros(M_.endo_nbr,M_.exo_nbr);
           oo_.dr.ghu(yp,:);
           zeros(M_.exo_nbr,M_.exo_nbr)];
    zxup = [zeros(M_.nspred,M_.nspred*M_.exo_nbr);
            zeros(M_.endo_nbr,M_.nspred*M_.exo_nbr);
            oo_.dr.ghxu(yp,:)*kron(oo_.dr.ghx(y_,:),Iu);
            zeros(M_.exo_nbr,M_.nspred*M_.exo_nbr)];
    zupup = [zeros(M_.nspred,M_.exo_nbr*M_.exo_nbr);
             zeros(M_.endo_nbr,M_.exo_nbr*M_.exo_nbr);
             oo_.dr.ghuu(yp,:);
             zeros(M_.exo_nbr,M_.exo_nbr*M_.exo_nbr)];
    zss = [zeros(M_.nspred,1);
           ghs2;
           oo_.dr.ghx(yp,:)*ghs2(y_,:) + ghs2(yp,:);
           zeros(M_.exo_nbr,1)];
    Fxupup = d1f_d1yp*oo_.dr.ghxuu(yp,:)*kron(oo_.dr.ghx(y_,:),Iu_kron_Iu) + d2f*(kron(zx,zupup) + kron(zup,zxup)*(P_up1_x1up2+P_up2_x1up1)) + d3f*kron(zx,kron(zup,zup));
    RHS = d1f_d1yp*oo_.dr.ghxx(yp,:)*kron(oo_.dr.ghx(y_,:),ghs2(y_,:)) + d2f*kron(zx,zss) + Fxupup*kron(Ix,SIGMA2);
    ghxss = gensylv(1,A,B,oo_.dr.ghx(y_,:),-RHS);
end

%% Re-compute ghuss
I_uuu = speye(M_.exo_nbr * M_.exo_nbr * M_.exo_nbr );
id_uuu = reshape(1:(M_.exo_nbr * M_.exo_nbr * M_.exo_nbr ),1,  M_.exo_nbr , M_.exo_nbr , M_.exo_nbr );
if ~isequal(SIGMA2,M_.Sigma_e(:))
    % permuted identity matrices useful for tensor unfolding, see e.g. Levintal (2017)    
    P_up1_u1up2 = I_uuu(:,ipermute(id_uuu,[0 flip(-[2,1,3]+4)]+1)); %this permutes e.g. kron(zup,zuup) to make it compatible with [zup]_d1[zuup]_b1d1 tensor
    P_up2_u1up1 = I_uuu(:,ipermute(id_uuu,[0 flip(-[3,1,2]+4)]+1)); %this permutes e.g. kron(zup,zuup) to make it compatible with [zup]_d2[zuup]_b1d2 tensor
    
    zu = [zeros(M_.nspred,M_.exo_nbr);
          oo_.dr.ghu;
          oo_.dr.ghx(yp,:)*oo_.dr.ghu(y_,:);
          Iu];
    zuup = [zeros(M_.nspred,M_.exo_nbr*M_.exo_nbr);
            zeros(M_.endo_nbr,M_.exo_nbr*M_.exo_nbr);
            oo_.dr.ghxu(yp,:)*kron(oo_.dr.ghu(y_,:),Iu);
            zeros(M_.exo_nbr,M_.exo_nbr*M_.exo_nbr)];
    Fuupup = d1f_d1yp*oo_.dr.ghxuu(yp,:)*kron(oo_.dr.ghu(y_,:),Iu_kron_Iu) + d2f*(kron(zu,zupup) + kron(zup,zuup)*(P_up1_u1up2+P_up2_u1up1)) + d3f*kron(zu,kron(zup,zup));
    RHS = d1f_d1yp*(oo_.dr.ghxx(yp,:)*kron(oo_.dr.ghu(y_,:),ghs2(y_,:)) + ghxss(yp,:)*oo_.dr.ghu(y_,:)) + d2f*kron(zu,zss) + Fuupup*kron(Iu,SIGMA2);
    ghuss = -A\RHS;
end

%% Re-compute ghs3
% permuted identity matrices useful for tensor unfolding, see e.g. Levintal (2017)
P_u1_u2u3 = I_uuu(:,ipermute(id_uuu,[0 flip(-[1,2,3]+4)]+1)); %this permutes e.g. kron(zu,zuu) to make it compatible with [zu]_b1[zuu]_b2b3 tensor
P_u2_u1u3 = I_uuu(:,ipermute(id_uuu,[0 flip(-[2,1,3]+4)]+1)); %this permutes e.g. kron(zu,zuu) to make it compatible with [zu]_b2[zuu]_b1b3 tensor
P_u3_u1u2 = I_uuu(:,ipermute(id_uuu,[0 flip(-[3,1,2]+4)]+1)); %this permutes e.g. kron(zu,zuu) to make it compatible with [zu]_b3[zuu]_b1b2 tensor
RHS = (d1f_d1yp*oo_.dr.ghuuu(yp,:) + d2f_d2yp*kron(oo_.dr.ghu(yp,:),oo_.dr.ghuu(yp,:))*(P_u1_u2u3+P_u2_u1u3+P_u3_u1u2) + d3f_d3yp*kron(oo_.dr.ghu(yp,:),kron(oo_.dr.ghu(yp,:),oo_.dr.ghu(yp,:))))*SIGMA3;
ghs3 = -((A+B)\RHS);

%% store into output structure
if ~isequal(SIGMA2,M_.Sigma_e(:))
    oo_.dr.ghs2 = ghs2;
    oo_.dr.ghxss = ghxss;
    oo_.dr.ghuss = ghuss;
end
oo_.dr.ghs3 = ghs3;
end % main function end
