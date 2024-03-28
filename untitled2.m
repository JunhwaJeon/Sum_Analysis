al=1/2; om=64;
fir_moment=SumMoment_4(al, om, 1);
sec_moment=SumMoment_4(al, om, 2);
four_moment=SumMoment_4(al, om, 4);
z1=log((fir_moment)^2/(sec_moment));
z2=log((sec_moment)^2/(four_moment));

% 이변수 함수 정의
syms m k
fun1=equa(m,k,z1);
fun2=equa(m,k,z2);
gradient(fun,[m,k])
fun3=(fun1)^2+fun2^2

function F=equa(m,k,z1)
F=g(m,k)-log(z1);
end

function result = g(m, k)
    h1 = 2*gammaln(m + 1/k);
    h2 = gammaln(m);
    h3 = gammaln(m + 2/k);
    
    result = h1 - h2 - h3;
end