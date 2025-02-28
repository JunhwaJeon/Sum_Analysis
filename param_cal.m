close all; clear; clc;

m=0:0.1:10000;
k=logspace(-3,-2.3,10^3);
al=1/2; om=64;
fir_moment=SumMoment_4(al, om, 1);
sec_moment=SumMoment_4(al, om, 2);
four_moment=SumMoment_4(al, om, 4);
z1=log((fir_moment)^2/(sec_moment));
z2=log((sec_moment)^2/(four_moment));

[M,K]=meshgrid(m,k);

Z1=((2*gammaln(M+1./K)-gammaln(M)-gammaln(M+2./K)-z1).^2);
Z2=((2*gammaln(M+2./K)-gammaln(M)-gammaln(M+4./K)-z2).^2);
sol=[0;0];


for i=1:length(m)
    for j=1:length(k)
        if (Z1(j,i)<=2333) && (Z2(j,i)<=2333)
            sol=[sol,[m(i);k(j)]];
        end
    end
end
sol;

mesh(M,K,Z1);