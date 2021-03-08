function [ul_list_sparse,old_ul_list,sMu,dlam,delul_list,free_dofs] = solve_arclength(it,delul_list,ul_list,old_ul_list,K,R,bc_dofs,ndofs,count,res,Fef,lambda,dlam,sMu)
%SOLVE_LINEAR_SYSTEM Summary of this function goes here
%   Detailed explanation goes here
phi=1.;
lp=0.5;
%% STORE OLD ONE
% old_ul_list=ul_list;
free_dofs=[1:ndofs]';
free_dofs(bc_dofs(:,1))=[];
ul_list_sparse=sparse(ndofs,1);

% ul_list_trial = zeros(size(ul_list));

%% METHOD 2: MAKE CORRESPONDING ROWS AND COLUMNS IDENTITY (KNOWN DISPLACEMNETS)
if count==1
    R=R-K(:,bc_dofs(:,1))*bc_dofs(:,2);
end
K(:,bc_dofs(:,1)) = [];
K(bc_dofs(:,1),:) = [];
% R(bc_dofs(:,1)) = [];
Fef(bc_dofs(:,1)) = [];
% res(bc_dofs(:,1)) = [];
%% Start Arclength
if count==1
    delres = zeros(size(free_dofs,1),1);
    delF   = K\Fef;
else
    delres = K\res;
    delF   = K\Fef;
end

%% compute root
c1 = dot(delF,delF) + phi.^2*dot(Fef,Fef);
c2 = 2*dot(delF,(delul_list+delres)) + 2*dlam*phi.^2*dot(Fef,Fef);
c3 = dot( (delul_list+delres) , (delul_list+delres) )-lp.^2+dlam.^2*phi.^2*dot(Fef,Fef);
DD = c2.^2 - 4*c1*c3;
if DD > 1e-13
    dellambda1 = (-c2+sqrt(DD))/(2*c1);
    dellambda2 = (-c2-sqrt(DD))/(2*c1);
    %     dellambda = max(dellambda1,dellambda2);
else
    error('no real root. Reduce the arclength. lp=lp/2 and reset the solution.')
end

%% determine direction

% %% predictor method 1
% if count==1
%     check = dot(old_ul_list(free_dofs),delF)+phi*dlam*dot(Fef,Fef);
%     % this is the first step
%     if check < 1e-14
%         if sign(dlam+dellambda1)==sign(eigs(K,1,'sm'))
%             dellambda=dellambda1;
%         else
%             dellambda=dellambda2;
%         end
%     else % end of first step
%         if sign(dlam+dellambda1)==sign(check)
%             dellambda=dellambda1;
%         else
%             dellambda=dellambda2;
%         end
%     end
% 
%     % now the corrector phase
% else
% 
%     %  check2= dot(dw,delF)+phi*phi*dlam;
%     %   minimum angle
% 
%     c1 = dot( delul_list , delul_list+delres+(dlam+dellambda1)*delF);
%     c2 = dot( delul_list , delul_list+delres+(dlam+dellambda2)*delF);
%     if c1 > c2
%         dellambda=dellambda1;
%     else
%         dellambda=dellambda2;
%     end
% 
% 
% end
%% predictor method 2


Att = eigs(K,5,'sm');

if it==1
    sMu=sign(prod(Att));
end

if sMu==sign(prod(Att)); %sign(eigs(Kt,1,'sm'))
    dellambda=dellambda1;
else
    dellambda=dellambda2;
end

%% method 3
% if count==1
%     check = dot(old_ul_list(free_dofs),delF)+phi*dlam*dot(Fef,Fef);
%     this is the first step
%     if check < 1e-14
%         if sign(dlam+dellambda1)==sign(eigs(K,1,'sm'))
%             dellambda=dellambda1;
%         else
%             dellambda=dellambda2;
%         end
%     else % end of first step
%         if sign(dlam+dellambda1)==sign(check)
%             dellambda=dellambda1;
%         else
%             dellambda=dellambda2;
%         end
%     end
% 
%     now the corrector phase
% else
% 
%     check2 = dot(delul_list,delF)+phi*phi*dlam;
%       minimum angle
%     if check2*dellambda1 > check2*dellambda2
%         dellambda=dellambda1;
%     else
%         dellambda=dellambda2;
%     end
% 
% 
% end



%% update solution
dlam = dlam + dellambda;
delw = delres+dellambda.*delF;
%
% update increment
delul_list = delul_list + delw;
% w(freeIndex) = w(freeIndex) + delw;
% w(uDof) = wold(uDof) + uFixed;
%     solved_system=K\R;
ul_list_sparse(free_dofs) = delw + ul_list(free_dofs);
ul_list_sparse(bc_dofs(:,1),:)=bc_dofs(:,2)+old_ul_list(bc_dofs(:,1),1);

% Ul is updated!
% compute new internal force
%
% fe =
%
% res = fe(free_dofs) + (lambda+dlam).*Fef;



end
