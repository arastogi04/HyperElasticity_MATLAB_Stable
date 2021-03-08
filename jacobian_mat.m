function [ jac_mat,det_jac_mat] = jacobian_mat( coord, dshape,debug,file01 )
%JACOBIAN_MAT Summary of this function goes here
%   Detailed explanation goes here

jac_mat=coord*dshape;
det_jac_mat=det(jac_mat);

if (det_jac_mat<0)
    fprintf('JACOBIAN NEGATIVE\n\n\n\nJACOBIAN NEGATIVE\n');
    fprintf(file01,'JACOBIAN NEGATIVE\n\n\n\nJACOBIAN NEGATIVE\n');
end


if(debug)
    mprint(jac_mat,'JACOBIAN MATRIX',file01);
    mprint(det_jac_mat,'DETERMINANT OF JACOBIAN MATRIX',file01);
end

end

