function     [d,sigma,sig,J_tilde, p_tilde, dPsi_vol_dJ_tilde, d2Psi_vol_dJ_tilde2,Je]=mate03d_vol_iso(F,J,mat_para,J_tilde,p_tilde,debug,file01)
% Summary of this function goes here
%   Detailed explanation goes here

% to store globally.

g12 = mat_para(1);
lambda = mat_para(2);
kappa= lambda + 2/3*g12;
Id2 = eye(3);
ttim = mat_para(4);
gp=mat_para(3);

g=1+gp*ttim;

Fe=F./g;
Je= det(Fe);
Jg = g*g*g;
%%

Nota= [ 1 4 6
        4 2 5
        6 5 3];

%%% J_tilde = v / V_0
%%% p = 1/V_0 \int tr(mixed sigma / (3 J_tilde))

be = Fe*Fe';
be_bar= Je^(-2/3)*be;

dPsi_vol_dJ_tilde= kappa/2*(J_tilde - 1/J_tilde);
d2Psi_vol_dJ_tilde2 = (kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde));
% fprintf('J_tilde %f Je_tilde %f \n',J_tilde, Je_tilde);
%%% TAU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TAU_VOL
taue_vol = (p_tilde*Je).*Id2;
%%% TAU_ISO
taue_bar = g12*be_bar;
taue_iso = taue_bar - 1/3*(taue_bar(1,1)+taue_bar(2,2)+taue_bar(3,3)).*Id2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr_taue_bar = taue_bar(1,1) +taue_bar(2,2) +taue_bar(3,3);

%%% Jc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                II4(i,j,k,l)=(0.5)*(Id2(i,k)*Id2(j,l) + Id2(i,l)*Id2(j,k));
                IIP(i,j,k,l) = II4(i,j,k,l) -1/3*(Id2(i,j)*Id2(k,l));
                %%% VOLUMETRIC PART OF TANGENT Jc_vol_tensor
                Jce_vol_tensor(i,j,k,l) = p_tilde *Je*(Id2(i,j)*Id2(k,l) - 2 *II4(i,j,k,l));
                
                %%% ISOCHORIC PART OF TANGENT Jc_iso_tensor
                Jce_iso_tensor(i,j,k,l) = 2/3*(tr_taue_bar)*IIP(i,j,k,l) - 2/3*(taue_iso(i,j)*Id2(k,l) +Id2(i,j)*taue_iso(k,l));
                
                Jce_vol(Nota(i,j),Nota(k,l))=Jce_vol_tensor(i,j,k,l);
                Jce_iso(Nota(i,j),Nota(k,l))=Jce_iso_tensor(i,j,k,l);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taue = taue_vol + taue_iso ;
Jce = Jce_vol + Jce_iso ;

sigma(1) = taue(1,1)/J;
sigma(2) = taue(2,2)/J;
sigma(3) = taue(3,3)/J;
sigma(4) = taue(1,2)/J;
sigma(5) = taue(2,3)/J;
sigma(6) = taue(1,3)/J;

sig=taue./J;
d = Jce./J ;

if(debug)
    mprint(sigma,'CAUCHY STRESS',file01);
    mprint(d,'CONSISTANT TANGENT',file01);
end

end

