function [forced_dofs,force_mag] = set_forc_boundary(mesh_input_file,define_BC,dof_list,fixed_dof)
%DISP_BOUNDARY_SET Summary of this function goes here
%   Detailed explanation goes here

%% OUTPUT VARIABLEs
% bc_dofs=[0 0];
forced_dofs=[];
force_mag=[];
%%
x_dof_cboun=[];
y_dof_cboun=[];
z_dof_cboun=[];
x_mag=[];
y_mag=[];
z_mag=[];

global ndim;
nodes=sprintf('input/%s.nodes',mesh_input_file);
data = dlmread(nodes);

if (isfield(define_BC,'CFORC'))
    
    CFORC=define_BC.CFORC; %=[2.0 -0.5 1 1];'
    
    for tt=1:size(CFORC,1)
        node_at_specified_coor(tt,1) = find(ismember(data(:,3:3+ndim-1),CFORC(tt,1:ndim),'rows'));
        
        if (CFORC(tt,ndim+1));
            x_dof_cboun(tt,:)=ndim*(node_at_specified_coor(tt,1)-1)+1;
            x_mag(tt,:) = CFORC(tt,ndim+1);
        end;
        if (CFORC(tt,ndim+2));
            y_dof_cboun(tt,:)=ndim*(node_at_specified_coor(tt,1)-1)+2;
            y_mag(tt,:) = CFORC(tt,ndim+2);
        end;
        if (ndim==3 && CFORC(tt,ndim+3));
            z_dof_cboun(tt,:)=ndim*(node_at_specified_coor(tt,1)-1)+3;
            z_mag(tt,:) = CFORC(tt,ndim+3);
        end;
    end
    
    chk4=isempty(x_dof_cboun);
    chk5=isempty(y_dof_cboun);
    if (ndim==3) chk6=isempty(z_dof_cboun); end;
    
    fixed_dofs_cboun=[]; % initilize
    fixed_mag_cboun=[]; % initilize

    if (~chk4)
        fixed_dofs_cboun=union(fixed_dofs_cboun,x_dof_cboun');
        fixed_mag_cboun = union(fixed_mag_cboun,x_mag');
    end
    if (~chk5)
        fixed_dofs_cboun=union(fixed_dofs_cboun,y_dof_cboun');
         fixed_mag_cboun = union(fixed_mag_cboun,y_mag');
    end
    if (ndim==3 && ~chk6)
        fixed_dofs_cboun=union(fixed_dofs_cboun,z_dof_cboun');
         fixed_mag_cboun = union(fixed_mag_cboun,z_mag');
    end
    
    fixed_nodes_cboun=node_at_specified_coor;
    
    forced_dofs=union(forced_dofs,fixed_dofs_cboun);
    force_mag=union(force_mag,fixed_mag_cboun);
    
end


% We have bc_dofs and fixed_dofs.
% OUTPUT
% bc_dofs=union(bc_dofs,fixed_dofs,'rows');


end

