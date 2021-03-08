function PlotMesh(coor_list_dof_or_coord,connect_list,fixed_nodes,disp_nodes)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Finite Element Method Mesh in current configuration
% Synopsis :
%           PlotMesh(coordinates,nodes,fixed_nodes,Disp_nodes)
%--------------------------------------------------------------------------
global ndim
global total_node_no

row_num=total_node_no;
col_num=size(coor_list_dof_or_coord,2);
coor_list=zeros(row_num,5);

if (ndim==3)
    
    if (col_num==1)
        coor_list(:,3)=coor_list_dof_or_coord(1:ndim:end,1);
        coor_list(:,4)=coor_list_dof_or_coord(2:ndim:end,1);
        coor_list(:,5)=coor_list_dof_or_coord(3:ndim:end,1);
    elseif (col_num==3)
        coor_list(:,3:5)=coor_list_dof_or_coord;
    elseif (col_num==5)
        coor_list=coor_list_dof_or_coord;
    end
    
    
    nel=size(connect_list,1);          %total element number
    total_node=size(coor_list,1);
    
    hold off;figure(1);hold on;view(3);set(0,'DefaultFigureVisible','off')
    
    for k=1:nel
        coordinates= coor_list(connect_list(k,4:end),:);
        
        % face vertices
        m(1,:)=[1 2 3 4];
        m(2,:)=[2 3 7 6];
        m(3,:)=[5 6 7 8];
        m(4,:)=[1 4 8 5];
        m(5,:)=[3 4 8 7];
        m(6,:)=[1 2 6 5];
        
        
        for t=1:6
            xyz=coordinates(m(t,:),3:end);
            fill3(xyz(:,1),xyz(:,2),xyz(:,3),connect_list(k,3),'FaceAlpha', 0.8);
        end
        
    end
    
    % plot DISP nodes.
    for i = 1:length(disp_nodes) % change to total node number
        disp_node=disp_nodes(i);
        quiver3(coor_list(disp_node,3),coor_list(disp_node,4),coor_list(disp_node,5),0.2,0.2,0.2)
    end
    % plot fixed nodes.
    for i = 1:length(fixed_nodes) % change to total node number
        fix_node=fixed_nodes(i);
        plot3(coor_list(fix_node,3),coor_list(fix_node,4),coor_list(fix_node,5),'.','markersize',20)
    end
    % write node number on figure;
    for i = 1:total_node % change to total node number
        text(coor_list(i,3),coor_list(i,4),coor_list(i,5), ['N' num2str(i)],'VerticalAlignment', 'top', 'FontSize', 12,'Color', 'r');
    end
    
    % TO DO: Quiver can be implemented for fixed notation.
    
    xlabel('x'); xlabel('y'); xlabel('z');
    axis equal;
    axis on;
    print('-dtiff', '-r200', ['output/plot_reference']);
    
elseif (ndim==2)
    
    
    if (col_num==1)
        coor_list(:,3)=coor_list_dof_or_coord(1:ndim:end,1);
        coor_list(:,4)=coor_list_dof_or_coord(2:ndim:end,1);
    elseif (col_num==2)
        coor_list(:,3:3+ndim-1)=coor_list_dof_or_coord;
    elseif (col_num==4)
        coor_list=coor_list_dof_or_coord;
    end
    
    
    nel=size(connect_list,1);          %total element number
    total_node=size(coor_list,1);
    
    hold off;figure(1);hold on;view(2);set(0,'DefaultFigureVisible','off')
    
    for k=1:nel
        coordinates= coor_list(connect_list(k,4:end),:);
        % face vertices
        m(1,:)=[1 2 3 4];
        xy=coordinates(m(1,:),3:end);
        fill(xy(:,1),xy(:,2),connect_list(k,3),'FaceAlpha', 0.8);
    end
    
    % plot fixed nodes.
    for i = 1:length(fixed_nodes) % change to total node number
        fix_node=fixed_nodes(i);
        plot(coor_list(fix_node,3),coor_list(fix_node,4),'.','markersize',20)
    end
    % plot DISP nodes.
    for i = 1:size(disp_nodes,1) % change to total node number
        disp_node=disp_nodes(i,1);
        norm_1=norm([disp_nodes(i,2) disp_nodes(i,3)]);
        disp_x=disp_nodes(i,2)/norm_1;
        disp_y=disp_nodes(i,3)/norm_1;
        quiver(coor_list(disp_node,3),coor_list(disp_node,4),disp_x,disp_y)
        plot(coor_list(disp_node,3),coor_list(disp_node,4),'.','color','r','markersize',20)
    end
    
    % write node number on figure;
    for i = 1:total_node % change to total node number
        text(coor_list(i,3),coor_list(i,4),['N' num2str(i)],'VerticalAlignment', 'top', 'FontSize', 12,'Color', 'r');
    end
    
    % TO DO: Quiver can be implemented for fixed notation.
    
    xlabel('x'); ylabel('y');
    axis equal;
    axis on;
    print('-dtiff', '-r200', ['output/plot_reference']);
    
    
end

end