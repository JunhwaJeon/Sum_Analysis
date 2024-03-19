function [alp_s, mu_s, om_s]=Weibull_MoM_4(al, om)
% Weibull을 따르는 Summand <4개> 의 <alpha>, <omega>를 입력받아
% Sum으로 Moment Matching하여 alpha, mu를 구하는 함수

fir_moment=SumMoment_4(al, om, 1);
sec_moment=SumMoment_4(al, om, 2);
four_moment=SumMoment_4(al, om, 4);


function F=equations(x,fir_moment,sec_moment, four_moment)
F(1)=(gamma(x(1)+1/x(2))^2)*(sec_moment/(fir_moment)^2)-gamma(x(1))*gamma(x(1)+2/x(2)); 
F(2)=(gamma(x(1)+2/x(2))^2)*(four_moment/(sec_moment)^2)-gamma(x(1))*gamma(x(1)+4/x(2));
end

% Solving with <fsolve>
fun=@(x)equations(x,fir_moment,sec_moment,four_moment);
options=optimoptions('fsolve','MaxFunctionEvaluations',400);
sol=fsolve(fun,[1 1], options);
alp_s=sol(1); mu_s=sol(2);
om_s=((mu_s^(1/alp_s))*gamma(mu_s)*fir_moment/(gamma(mu_s+1/alp_s)))^alp_s;


%{
% Solving with <vpasolve>
syms x1 x2
F=[(gamma(x1+1/x2)^2)*(sec_moment/(fir_moment)^2)==gamma(x1)*gamma(x1+2/x2)
,(gamma(x1+2/x2)^2)*(four_moment/(sec_moment)^2)==gamma(x1)*gamma(x1+4/x2)];
sol=vpasolve(F, [x1 x2]);
[alp_s, mu_s]=sol(:);
om_s=((mu_s^(1/alp_s))*gamma(mu_s)*fir_moment/(gamma(mu_s+1/alp_s)))^alp_s;
%}
end