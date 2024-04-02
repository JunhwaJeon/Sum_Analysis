cvx_clear
cvx_solver('sedumi')
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
   minimize(2 * ((m + k - 1/2) * log(m + k) + 1/(12*(m + k)) - 1/(360*(m + k)^3) + 1/(1260*(m + k)^5) - 1/(1680*(m + k)^7)) ...
        - (m - 1/2) * log(m) + 1/(12*m) - 1/(360*m^3) + 1/(1260*m^5) - 1/(1680*m^7) ...
        - ((m + 2*k) - 1/2) * log(m + 2*k) + 1/(12*(m + 2*k)) - 1/(360*(m + 2*k)^3) + 1/(1260*(m + 2*k)^5) - 1/(1680*(m + 2*k)^7) ...
        - z1)
    %{
    minimize((2*((m + k-1/2)*log(m + k)+1/(12*(m + k))-1/(360*(m + k)^3)+1/(1260*(m + k)^5)-1/(1680*(m + k)^7))- ...
        (m-1/2)*log(m)+1/(12*m)-1/(360*m^3)+1/(1260*m^5)-1/(1680*m^7)- ...
        ((m + 2*k)-1/2)*log(m + 2*k)+1/(12*(m + 2*k))-1/(360*(m + 2*k)^3)+1/(1260*(m + 2*k)^5)-1/(1680*(m + 2*k)^7)-z1)^2+ ...
        (2*((m + k/2-1/2)*log(m + k/2)+1/(12*(m + k/2))-1/(360*(m + k/2)^3)+1/(1260*(m + k/2)^5)-1/(1680*(m + k/2)^7))- ...
        (m-1/2)*log(m)+1/(12*m)-1/(360*m^3)+1/(1260*m^5)-1/(1680*m^7)- ...
        ((m + k)-1/2)*log(m + k)+1/(12*(m + k))-1/(360*(m + k)^3)+1/(1260*(m + k)^5)-1/(1680*(m + k)^7)-z2)^2)
    %}
    subject to
        m > 0;
        k > 0;
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