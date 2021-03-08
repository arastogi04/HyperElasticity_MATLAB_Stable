function [K,Fd] = elmt13_par(d,ul,xl,node_p_el,int_rule,total_node_no,ndim,element_type,node_p_elm,debug,file01);
% MAIN PURPOSE: GET STIFFNESS AND RESIDUAL VALUE (s,r)

% xl
% ul
%     Get quadrature information - EVALUATE POSITIONS AND W8S
[int_points,int_weights]=gauss_int(int_rule,debug,file01);


%     Compute current geometry (ndimXnnode)
x_reff=xl;
x_curr=x_reff+ul;

%  ---------------------------------------------------------------------
%
%      LOOP OVER THE INTGRATION POINTS
%
%  ---------------------------------------------------------------------
K_mat=zeros(node_p_elm*ndim,node_p_elm*ndim);
K_geo_exp=zeros(node_p_elm*ndim,node_p_elm*ndim);
K_geo_pseudo=zeros(node_p_elm,node_p_elm);
Fd=zeros(node_p_elm*ndim,1);
for intip=1:node_p_el
    % Evaluate shape functions & derivatives
    [N_shape, dNdxi_shape]=shape_fun(int_points(:,intip),element_type,debug,file01);
    
    % Compute Jacobian: REFERENCE CONFIGURATION
    [jac_ref,det_jac_ref]=jacobian_mat(x_reff',dNdxi_shape,debug,file01);
    
    % Compute Jacobian: CURRENT CONFIGURATION
    [jac_cur,det_jac_cur]=jacobian_mat(x_curr',dNdxi_shape,debug,file01);
    
    % Compute deformation gradient F = je*Je^-1
    F=jac_ref\jac_cur;
    det_F=det_jac_cur/det_jac_ref;
    
    % Compute B-operator
    %       B-operator looks like that
    %       | Nk,x    0       0|
    %       |   0    Nk,y     0|
    %       |   0     0    Nk,z|
    %       | Nk,y   Nk,x     0|
    %       |    0   Nk,z  Nk,y|
    %       | Nk,z    0    Nk,x|
    
    inv_jac_cur=inv(jac_cur);
    dNidx=inv_jac_cur'*dNdxi_shape';
    
    if (ndim==3)
        B=zeros(6,3*node_p_elm);
        B(1,1:ndim:ndim*node_p_elm)=dNidx(1,:);
        B(2,2:ndim:ndim*node_p_elm)=dNidx(2,:);
        B(3,3:ndim:ndim*node_p_elm)=dNidx(3,:);
        
        B(4,1:ndim:ndim*node_p_elm)=dNidx(2,:);
        B(4,2:ndim:ndim*node_p_elm)=dNidx(1,:);
        
        B(5,2:ndim:ndim*node_p_elm)=dNidx(3,:);
        B(5,3:ndim:ndim*node_p_elm)=dNidx(2,:);
        
        B(6,1:ndim:ndim*node_p_elm)=dNidx(3,:);
        B(6,3:ndim:ndim*node_p_elm)=dNidx(1,:);
        
    elseif (ndim==2)
        B=zeros(3,2*node_p_elm);
        
        B(1,1:ndim:ndim*node_p_elm)=dNidx(1,:);
        B(2,2:ndim:ndim*node_p_elm)=dNidx(2,:);
        
        B(3,1:ndim:ndim*node_p_elm)=dNidx(2,:);
        B(3,2:ndim:ndim*node_p_elm)=dNidx(1,:);
    end
    
    % compute infinitesimal volume element
    dv_cur = det_jac_cur*int_weights(intip);
    %%% MATERIAL MODEL
    %     [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01(F,det_F,d,debug,file01);
    %     [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01_growing(F,det_F,d,debug,file01);
    [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01_2d(F,det_F,d,debug,file01);
    
    % Compute MATerial STIFFNESS matrix
    K_mat=K_mat+B'*cons_tangent_voigt*B*dv_cur;
    
    % Compute RESIDUAL
    Fd=Fd-B'*cauchy_sigma_voigt'*dv_cur;
    
    % Compute GEOmetric STIFFNESS matrix
    K_geo_pseudo=K_geo_pseudo+dNidx'*cauchy_sigma*dNidx*dv_cur;
end

if (ndim==3)
    for tt=0:node_p_elm-1
        K_geo_exp(1+3*tt,1:3:3*node_p_elm)=K_geo_pseudo(tt+1,:);
        K_geo_exp(2+3*tt,2:3:3*node_p_elm)=K_geo_pseudo(tt+1,:);
        K_geo_exp(3+3*tt,3:3:3*node_p_elm)=K_geo_pseudo(tt+1,:);
    end
elseif (ndim==2)
    for tt=0:node_p_elm-1
        K_geo_exp(1+2*tt,1:2:2*node_p_elm)=K_geo_pseudo(tt+1,:);
        K_geo_exp(2+2*tt,2:2:2*node_p_elm)=K_geo_pseudo(tt+1,:);
    end
end


K=K_mat+K_geo_exp;

end