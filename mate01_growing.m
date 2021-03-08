function     [aa,sigmae,sige]=mate01_growing(F,J,d,debug,file01)
%MATE01 Summary of this function goes here
%   Detailed explanation goes here

g12 = d(1,1);
lambda = d(2,1);
gp=d(3,1);
ttim=d(4,1);

% ISOTROPIC GROWTH
g=1+gp*ttim;

iso_flg=1;
ani_flg=0;

if (iso_flg)
    Fe=F./g;
elseif (ani_flg)
    n_0=[0;1;0];
    n_1=F*n_0;
    ontheta=1/g;
    Iani=n_1*n_0';
    Fe=F./ontheta +(1 - ontheta).*Iani;
else
    disp('ISOTROPIC OR ANISOTROPIC');
end

Je=det(Fe);
be = Fe*Fe';

kirche = g12.*be+(-g12+lambda/2*Je*Je-lambda/2).*eye(3);

sigmae(1) = kirche(1,1)/J;
sigmae(2) = kirche(2,2)/J;
sigmae(3) = kirche(3,3)/J;
sigmae(4) = kirche(1,2)/J;
sigmae(5) = kirche(2,3)/J;
sigmae(6) = kirche(1,3)/J;

sige=kirche./J;

Id2 = [1 0 0;0 1 0; 0 0 1];

Nota= [ 1 4 6
    4 2 5
    6 5 3 ];

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                II4(i,j,k,l)=(0.5)*(Id2(i,k)*Id2(j,l) + Id2(i,l)*Id2(j,k));
                CmatE(i,j,k,l)= lambda.*Je*Je*Id2(i,j)*Id2(k,l)/J + (2*g12 + lambda - lambda*Je*Je)*II4(i,j,k,l)/J;
                aa(Nota(i,j),Nota(k,l))=CmatE(i,j,k,l);
            end
        end
    end
end

if(debug)
    mprint(sigmae,'CAUCHY STRESS',file01);
    mprint(aa,'CONSISTANT TANGENT',file01);
end


end

