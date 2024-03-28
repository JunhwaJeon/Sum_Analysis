clear all; close all; clc;

al=1/2; om=64;
fir_moment=SumMoment_4(al, om, 1);
sec_moment=SumMoment_4(al, om, 2);
four_moment=SumMoment_4(al, om, 4);
z1=log((fir_moment)^2/(sec_moment));
z2=log((sec_moment)^2/(four_moment));
syms m k
objective=(g(m,k)-z1)^2+(g(m,k/2)-z2)^2;

cvx_begin
    variables m k
    minimize(objective(m,k));
    subject to
        m>0;
        k>0;
cvx_end

disp(['최적해: x = ', num2str(m), ', y = ', num2str(k)]);
disp(['최적값: ', num2str(cvx_optval)]);

function F=equa(m,k,z1,z2)
F=(g(m,k)-z1)^2+(g(m,k/2)-z2)^2;
end

function result = g(m, k)
    h1 = 2*f(m + 1/k);
    h2 = f(m);
    h3 = f(m + 2/k);
    
    result = h1 - h2 - h3;
end

function result=f(x)
    result=(x-1/2)*log(x)+1/(12*x)-1/(360*x^3)+1/(1260*x^5)-1/(1680*x^7);
end