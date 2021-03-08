function [fixed_dofs,fixed_nodes] = set_dirichlet_boundary(mesh_input_file,define_BC,dof)

global ndim

%% OUTPUT VALUES
fixed_dofs=[]; % initilize
fixed_nodes=[];
%%

nodes=sprintf('input/%s.nodes',mesh_input_file);

data = dlmread(nodes);

%% EBOUN implementation
if (isfield(define_BC,'EBOUN'))
    
    EBOUN=define_BC.EBOUN; %=[2 -0.5 1 1];
    
    for tt=1:size(EBOUN,1)
        nodes_at_edge1(tt,:)=find(data(:,2+EBOUN(tt,1))==EBOUN(tt,2));
        if (EBOUN(tt,3)); x_dof(tt,:)=ndim*(nodes_at_edge1(tt,:)-1)+1;  end;
        if (EBOUN(tt,4)); y_dof(tt,:)=ndim*(nodes_at_edge1(tt,:)-1)+2;  end;
        if (ndim==3 && EBOUN(tt,5)); z_dof(tt,:)=ndim*(nodes_at_edge1(tt,:)-1)+3;  end;
    end
    
    chk1=isempty(x_dof);
    chk2=isempty(y_dof);
    if (ndim==3) chk3=isempty(z_dof); end;
    
    %     fixed_dofs=[]; % initilize
    
    if (~chk1)
        fixed_dofs=union(fixed_dofs,x_dof');
    end
    if (~chk2)
        fixed_dofs=union(fixed_dofs,y_dof');
    end
    if (ndim==3 && ~chk3)
        fixed_dofs=union(fixed_dofs,z_dof');
    end
    
    fixed_nodes=sort(reshape(nodes_at_edge1,[],1));
end
%% CBOUN implementation
if (isfield(define_BC,'CBOUN')) 
    
    CBOUN=define_BC.CBOUN; %=[2.0 -0.5 1 1];'
    
    for tt=1:size(CBOUN,1)
        node_at_specified_coor(tt,1) = find(ismember(data(:,3:3+ndim-1),CBOUN(tt,1:ndim),'rows'));
        if (CBOUN(tt,ndim+1)); x_dof_cboun(tt,:)=ndim*(node_at_specified_coor(tt,1)-1)+1;  end;
        if (CBOUN(tt,ndim+2)); y_dof_cboun(tt,:)=ndim*(node_at_specified_coor(tt,1)-1)+2;  end;
        if (ndim==3 && CBOUN(tt,ndim+3)); z_dof_cboun(tt,:)=ndim*(node_at_specified_coor(tt,1)-1)+3;  end;
    end
    
    chk4=isempty(x_dof_cboun);
    chk5=isempty(y_dof_cboun);
    if (ndim==3) chk6=isempty(z_dof_cboun); end;
    
    fixed_dofs_cboun=[]; % initilize
    
    if (~chk4)
        fixed_dofs_cboun=union(fixed_dofs_cboun,x_dof_cboun');
    end
    if (~chk5)
        fixed_dofs_cboun=union(fixed_dofs_cboun,y_dof_cboun');
    end
    if (ndim==3 && ~chk6)
        fixed_dofs_cboun=union(fixed_dofs_cboun,z_dof_cboun');
    end
    
    fixed_nodes_cboun=node_at_specified_coor;
    
    fixed_dofs=union(fixed_dofs,fixed_dofs_cboun);
    fixed_nodes=union(fixed_nodes,fixed_nodes_cboun);
    
end

end
