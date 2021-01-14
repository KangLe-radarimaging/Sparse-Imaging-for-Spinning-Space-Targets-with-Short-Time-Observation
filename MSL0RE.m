function ppp=MSL0RE(scc,T,T_inv)
% scc:             signal to be reconstructed
% T=Phi*Psi:       measurement matrix
% T_inv:           pseudo inverse of the measurement matrix
[~,N]=size(T);
nn1=length(scc);
y=(T.').^(-1)*scc;
yy=2*max(abs(y));
aug_y=y;
while yy>0.01
    r_n=aug_y;
yyy=yy;
for l=1:10
    theta=r_n.*exp(-abs(r_n).^2/(yyy^2));
    r_n=r_n-2*theta;
    r_n=r_n-T_inv*(T*r_n-scc);
end
yy=0.5*yy;
aug_y=r_n;
end
ppp=aug_y;

end
