function [K,Fd,cell_stored] = mixed_elmt13(d,ul,xl,cell_stored)
% MAIN PURPOSE: GET STIFFNESS AND RESIDUAL VALUE (s,r)
% s=zeros(total_node_no,total_node_no); r=zeros(total_node_no,1);
% mat_para=[8,392,0]; % [g12,lambda,gc]
% mat_para=[8,392,0]; % [g12,lambda,gc]

% _________________ ELMT13.m ____________________
global node_p_el int_rule debug file01 total_node_no ndim element_type node_p_elm

%     Get quadrature information - EVALUATE POSITIONS AND W8S
if (debug); fprintf(file01,'\n GAUSS QUADRATURE - elmt13.m \n '); end;
[int_points,int_weights]=gauss_int(int_rule,debug,file01);
% ul
% xl

%     Compute current geometry (ndimXnnode)
x_reff=xl;
x_curr=x_reff+ul;

if (debug)
    mprint(x_reff','REFERENCE COORDINATES (X,Y,Z)',file01);
    mprint(ul','Deformation vector (X,Y,Z)',file01);
    mprint(x_curr','SPATIAL COORDINATES (X,Y,Z)',file01);
end

%  ---------------------------------------------------------------------
%
%      LOOP OVER THE INTGRATION POINTS
%
%  ---------------------------------------------------------------------
K_mat=zeros(node_p_elm*ndim,node_p_elm*ndim);
K_geo=zeros(node_p_elm*ndim,node_p_elm*ndim);
K_geo_reduced=zeros(node_p_elm,node_p_elm);
Fd=zeros(node_p_elm*ndim,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_up = zeros(24,1);
k_Jp = zeros(1,1);
k_JJ = zeros(1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% SHAPE FUNCTIONS FOR J AND p FIELD.
N_J=1;
N_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for intip=1:size(int_points,2)
    if (debug); fprintf(file01,'\n (%d)^th LOOP OVER THE INTGRATION POINTS - elmt13.m  _____ \n ',intip); end;
    
    % Evaluate shape functions & derivatives
    [N_shape, dNdxi_shape]=shape_fun(int_points(:,intip),element_type,debug,file01);
    
    if (debug); fprintf(file01,'\nREFERENCE'); end;
    % Compute Jacobian: REFERENCE CONFIGURATION
    [jac_ref,det_jac_ref]=jacobian_mat(x_reff',dNdxi_shape,debug,file01);
    
    if (debug); fprintf(file01,'\nCURRENT'); end;
    % Compute Jacobian: CURRENT CONFIGURATION
    [jac_cur,det_jac_cur]=jacobian_mat(x_curr',dNdxi_shape,debug,file01);
    
    % Compute deformation gradient F = je*Je^-1
    F=jac_cur*inv(jac_ref);
    %     F=jac_ref\jac_cur; ! NOT THE SAME THING
    det_F=det_jac_cur/det_jac_ref;
    
    if (debug)
        mprint(jac_ref,'jac_ref ',file01);
        mprint(jac_cur,'jac_cur ',file01);
        mprint(F,'DEFORMATION GRADIENT',file01);
        mprint(det_F,'DETERMINANT OF DEFORMATION GRADIENT',file01);
        mprint(det_jac_ref,'det_jac_ref ',file01);
        mprint(det_jac_cur,'det_jac_cur ',file01);
    end
    
    if(debug)
        mprint(F,'DEFORMATION GRADIENT',file01);
        mprint(det_F,'DETERMINANT OF DEFORMATION GRADIENT',file01);
    end
    
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
    
    if (debug)
        mprint(B,'B OPERATOR',file01);
        mprint(dNidx,'DERIVATIVES OF SHAPE FUNCTIONS WITH RESPECT TO CARTESIAN COORDINATES',file01);
    end
    
    % compute infinitesimal volume element
    dv_cur = det_jac_cur*int_weights(intip);
    dv_ref = det_jac_ref*int_weights(intip);
    %%% MATERIAL MODEL
    
    %         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01(F,det_F,d,debug,file01);
    %         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01_2d(F,det_F,d,debug,file01);
    %         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01_growing(F,det_F,d,debug,file01);
    [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma,cell_stored]=mate03d_vol_iso(F,det_F,d,debug,file01,cell_stored);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% COMPUTE K_uu (24X24) 

    % Compute MATerial STIFFNESS matrix
    K_mat=K_mat+B'*cons_tangent_voigt*B*dv_cur;
    % Compute GEOmetric STIFFNESS matrix
    K_geo_reduced=K_geo_reduced+dNidx'*cauchy_sigma*dNidx*dv_cur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% COMPUTE K_pu (24X1) -> when pressure has the shape function -1-. 
    m = [1 1 1 0 0 0]';
    k_up=k_up+B'*m*N_p*dv_cur;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% COMPUTE K_jp (1X1) -> when pressure has the shape function -1-. 
    k_Jp=k_Jp+N_J'*N_p*dv_ref;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% COMPUTE K_jj (1X1) -> when pressure has the shape function -1-. 
    d2Psi_vol_dJ_tilde2 = cell_stored{3,2}; 
    
    k_JJ = k_JJ + N_J'* d2Psi_vol_dJ_tilde2 * N_J* dv_ref;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute RESIDUAL
    Fd=Fd-B'*cauchy_sigma_voigt'*dv_cur;
end

if (ndim==3)
    for tt=0:node_p_elm-1
        K_geo(1+3*tt,1:3:3*node_p_elm)=K_geo_reduced(tt+1,:);
        K_geo(2+3*tt,2:3:3*node_p_elm)=K_geo_reduced(tt+1,:);
        K_geo(3+3*tt,3:3:3*node_p_elm)=K_geo_reduced(tt+1,:);
    end
elseif (ndim==2)
    for tt=0:node_p_elm-1
        K_geo(1+2*tt,1:2:2*node_p_elm)=K_geo_reduced(tt+1,:);
        K_geo(2+2*tt,2:2:2*node_p_elm)=K_geo_reduced(tt+1,:);
    end
end

K=K_mat+K_geo;



end