clear all; close all; clc;

al=1/2; om=64;
fir_moment=SumMoment_4(al, om, 1);
sec_moment=SumMoment_4(al, om, 2);
four_moment=SumMoment_4(al, om, 4);
z1=log((fir_moment)^2/(sec_moment));
z2=log((sec_moment)^2/(four_moment));

x=optimvar('x',2,'LowerBound',[0,0.1^10]);
eq1=2*gammaln(x(1) + 1/x(2))-gammaln(x(1) + 2/x(2))-gammaln(x(1))-z1==0;
eq2=2*gammaln(x(1) + 2/x(2))-gammaln(x(1) + 4/x(2))-gammaln(x(1))-z2==0;

prob=eqnproblem;
prob.Equations.eq1=eq1;
prob.Equations.eq2=eq2;

x0.x=[0,1];
[sol,fval,exitflag]=solve(prob,x0);

function result = g(m, k)
    h1 = 2*gammaln(m + 1/k);
    h2 = gammaln(m);
    h3 = gammaln(m + 2/k);
    
    result = h1 - h2 - h3;
end