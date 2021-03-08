function [ul_list_sparse,old_ul_list,R] = solve_LS(ul_list,old_ul_list,K,R,bc_dofs,ndofs,count)
%SOLVE_LINEAR_SYSTEM Summary of this function goes here
%   Detailed explanation goes here

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
R(bc_dofs(:,1)) = [];

solved_system=K\R;

ul_list_sparse(free_dofs) = solved_system + ul_list(free_dofs);
ul_list_sparse(bc_dofs(:,1),:)=bc_dofs(:,2)+old_ul_list(bc_dofs(:,1),1);


end
