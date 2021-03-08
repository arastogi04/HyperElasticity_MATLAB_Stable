function [int_points,int_weights]=gauss_int(int_rule,debug,file01)
global ndim

if (strcmp(int_rule,'3x3'))
    
    p1 = sqrt(0.6);
    int_points=p1.*[-1  1  1 -1  0  1  0 -1  0
        -1 -1  1  1 -1  0  1  0  0 ];
    
    int_weights=[ 25/81 25/81 25/81 25/81 40/81 40/81 40/81 40/81 64/81];
    
elseif (strcmp(int_rule,'2x2'))
    
    p1 = sqrt(1/3);
    int_points=p1.*[-1  1  1 -1
        -1 -1  1  1 ];
    
    int_weights=[1 1 1 1];
    
elseif (strcmp(int_rule,'2x2x2'))
    p1 = 1/sqrt(3);
    int_points=p1.*[-1  1 -1  1 -1  1 -1  1
        -1 -1  1  1 -1 -1  1  1
        -1 -1 -1 -1  1  1  1  1];
    int_weights=zeros(1,8)+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(int_rule,'1point'))
    int_points=[0
        0];
    if (ndim==3) int_points(3,1)=0; end
    int_weights=4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif (strcmp(int_rule,'3x3x3'))
    
    
    x1D(1) = -sqrt(3.0/5.0);
    x1D(2) = 0.;
    x1D(3) = sqrt(3.0/5.0);
    
    % write(*,'(3f22.17 /)') (x1D(i), i = 1:3)
    % x1D = [-0.7745966692d0,zero,0.7745966692d0];
    for k = 1:3
        for j = 1:3
            for i = 1:3
                n = 9*(k-1) + 3*(j-1) + i;
                int_points(1,n) = x1D(i);
                int_points(2,n) = x1D(j);
                int_points(3,n) = x1D(k);
            end
        end
    end
    
    w1D(1) = 5.d0/9.d0;
    w1D(2) = 8.d0/9.d0;
    w1D(3) = 5.d0/9.d0 ;
    
    %!w1D = [0.555555555d0,0.888888888d0,0.55555555555d0]
    for k = 1:3
        for j = 1:3
            for i = 1:3
                n = 9*(k-1)+3*(j-1)+i;
                int_weights(n) = w1D(i)*w1D(j)*w1D(k);
            end
        end
    end
    
    
    
    
else
    if(debug); fprintf(file01,'ERROR in elmt13.f - wrong quadrature rule'); end;
end

if(debug)
    fprintf(file01,'\nQuadrature rule 3x3 \n');
    mprint(int_points,'LOCATIONS OF INTEGRATION POINTS',file01);
    mprint(int_weights,'WEIGTHS OF INTEGRATION POINTS',file01);
end


end