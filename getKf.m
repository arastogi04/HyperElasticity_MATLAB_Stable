function [K,R] = getKf(d,ul,xl,connect_list,add_field_e,element_dof,ndofs,sparse_sum)
global total_node_no
global ndim
global total_element_no

K = sparse(total_node_no*ndim,total_node_no*ndim);
R = sparse(total_node_no*ndim,1);
%% LOOP ALL OVER ELEMENTS
% SERIAL
for ie=1:total_element_no
    
    % ELEMENT FORMULATION WILL COME HERE.
    % INPUTS: ul: deformation vector
    %       : xl: reference coordinates
    
    [Ke,fe]=elmt13(d,ul(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1));
    %               [Ke,fe,add_field_e]=mixed_elmt13(d,ul(connect_list(ie,4:end)',1:ndim),du(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1),add_field(:,ie));
    %                 [Ke,fe,add_field_e]=mixed_elmt13_no_growth(d,ul(connect_list(ie,4:end)',1:ndim),du(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1),add_field(:,ie));
    
    add_field(:,ie)=add_field_e;
    %% ASSEMBLY PART
    [row_map col_map index_dofs] = node_mapping_matrix(element_dof(ie,:));
    K = sparse_sum(K,sparse(row_map,col_map,Ke,ndofs,ndofs));
    R = sparse_sum(R,sparse(index_dofs,1,fe,ndofs,1));
end
end