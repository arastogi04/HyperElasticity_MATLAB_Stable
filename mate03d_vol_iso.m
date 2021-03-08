function     [d,sigma,sig,J_tilde, p_tilde, dPsi_vol_dJ_tilde, d2Psi_vol_dJ_tilde2]=mate03d_vol_iso(F,J,mat_para,J_tilde,p_tilde,debug,file01)
% Summary of this function goes here
%   Detailed explanation goes here

% to store globally.

g12 = mat_para(1);
lambda = mat_para(2);
kappa= lambda + 2/3*g12;
Id2 = eye(3);

Nota= [ 1 4 6
    4 2 5
    6 5 3];

%%% J_tilde = v / V_0
%%% p = 1/V_0 \int tr(mixed sigma / (3 J_tilde))

b = F*F';
b_bar= J^(-2/3)*b;

dPsi_vol_dJ_tilde= kappa/2*(J_tilde - 1/J_tilde);
d2Psi_vol_dJ_tilde2 = (kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde));

%%% TAU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TAU_VOL
tau_vol = (p_tilde*J).*Id2;
%%% TAU_ISO
tau_bar = g12*b_bar;
tau_iso = tau_bar - 1/3*(tau_bar(1,1)+tau_bar(2,2)+tau_bar(3,3)).*Id2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr_tau_bar = tau_bar(1,1) +tau_bar(2,2) +tau_bar(3,3);

%%% Jc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                II4(i,j,k,l)=(0.5)*(Id2(i,k)*Id2(j,l) + Id2(i,l)*Id2(j,k));
                IIP(i,j,k,l) = II4(i,j,k,l) -1/3*(Id2(i,j)*Id2(k,l));
                %%% VOLUMETRIC PART OF TANGENT Jc_vol_tensor
                Jc_vol_tensor(i,j,k,l) = p_tilde * J*(Id2(i,j)*Id2(k,l) - 2 *II4(i,j,k,l));
                
                %%% ISOCHORIC PART OF TANGENT Jc_iso_tensor
                Jc_iso_tensor(i,j,k,l) = 2/3*(tr_tau_bar)*IIP(i,j,k,l) - 2/3*(tau_iso(i,j)*Id2(k,l) +Id2(i,j)*tau_iso(k,l));
                
                Jc_vol(Nota(i,j),Nota(k,l))=Jc_vol_tensor(i,j,k,l);
                Jc_iso(Nota(i,j),Nota(k,l))=Jc_iso_tensor(i,j,k,l);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = tau_vol + tau_iso ;
Jc = Jc_vol + Jc_iso ;

sigma(1) = tau(1,1)/J;
sigma(2) = tau(2,2)/J;
sigma(3) = tau(3,3)/J;
sigma(4) = tau(1,2)/J;
sigma(5) = tau(2,3)/J;
sigma(6) = tau(1,3)/J;

sig=tau./J;
d = Jc./J ;

if(debug)
    mprint(sigma,'CAUCHY STRESS',file01);
    mprint(d,'CONSISTANT TANGENT',file01);
end

end

