function [coor_list,connect_list,element_dof,dof] = Mesh_main(mesh_input_file)
%
%
%  (MESH + BOUNDARY C. + LOADING) DESCRIPTIONS + GRAPHICAL OUTPUT
%
%  1) GET THE MESH AS FEAP GET FROM CUBIT
global debug total_element_no total_node_no file01 node_p_elm ndim element_type
nodes=sprintf('input/%s.nodes',mesh_input_file);
elems=sprintf('input/%s.elems',mesh_input_file);

% FOR 3D MESH, READ NODE COORDINATES
% [node_no,dummy,x_coor,y_coor,z_coor] = coor_list
coor_list = dlmread(nodes); % get this from main program. main_FES.m
connect_list = dlmread(elems);

total_element_no=size(connect_list,1);          %total element number
total_node_no=size(coor_list,1);
node_p_elm=size(connect_list,2)-3;

dof=zeros(total_node_no*ndim,1);
element_dof=zeros(total_element_no,node_p_elm*ndim);

if (ndim==3)
    
    x=ndim*(connect_list(:,4:end)-1)+1;
    y=ndim*(connect_list(:,4:end)-1)+2;
    z=ndim*(connect_list(:,4:end)-1)+3;
    
    element_dof(:,1:3:end)=x(:,:);
    element_dof(:,2:3:end)=y(:,:);
    element_dof(:,3:3:end)=z(:,:);
    
    dof(1:ndim:end,1)=coor_list(:,3);
    dof(2:ndim:end,1)=coor_list(:,4);
    dof(3:ndim:end,1)=coor_list(:,5);
    
elseif (ndim==2)
        
    x=ndim*(connect_list(:,4:end)-1)+1;
    y=ndim*(connect_list(:,4:end)-1)+2;
    
    element_dof(:,1:2:end)=x(:,:);
    element_dof(:,2:2:end)=y(:,:);
    
    dof(1:ndim:end,1)=coor_list(:,3);
    dof(2:ndim:end,1)=coor_list(:,4);
    
end

% coor_list_check(:,1)=dof(1:ndim:end,1)
% coor_list_check(:,2)=dof(2:ndim:end,1)
% coor_list_check(:,3)=dof(3:ndim:end,1)

% COMPATIBILITY FOR FEAP
% coor_list = [NODE,0, X   Y   Z]
% connect_list = [ELEMENT,0,MAT_NO,CONNECTIVITY OF NODES]
fprintf(file01,'\nElement type: %10s \t Node Number per element: %2d \nTotal node number: %6d \t Total element number: %6d\n\n ',element_type,node_p_elm,total_node_no,total_element_no);


end