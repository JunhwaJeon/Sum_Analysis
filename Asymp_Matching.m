clear all;
close all;
clc;

%% Channel Matrix Under Spatial Correlation Fromulation%%
sq2 = sqrt(0.5); %상수 지정

%iid Rayleigh channel matrix
nR=8; nT=8; %Antenna 개수 지정
H_w=sq2*(randn(nR,nT)+1j*randn(nR,nT));

%spatial correlation matrix
%random correlation coefficient 주어서 어떤 distribution이 matching이 잘 되는지 확인
phi_T=zeros(nT);
for i=1:nT
    phi_T(i,i)=1;
    for j=i+1:nT
        x=-1+2*rand;
        phi_T(i,j)=x+(2*randi(2)-3)*1j*sqrt(1-x^2);
        phi_T(j,i)=conj(phi_T(i,j));
    end
end
t_corr=sqrtm(phi_T); %transmit correlation mtx
H=H_w*t_corr; %channel matrix

%% Define Exact Sum %%

%scale parameter 구하기
sum=0;
for i=1:nT
    sum=sum+t_corr(i,2);
end
scl_para=abs(sum)^4;

%Weibull distribution of shape parameter=1/2 정의하기
%pdf_wei(x)=(1/(scl_para*2))*((x/scl_para)^(-1/2))*exp(-(x/scl_para)^(1/2));
%Weibull sum 정의하기 - 다중적분 이용하여 정의
%sum_wei(x)=

%% Define Original Approximation %%

%% Define Asymptotic Matching Approximation %%

%summand's x=0 Taylor series expansion
%Summand's coefficient are given by: a=1/(2sqrt(t)), b=-1/2
%Sum coefficient are given by: a=((sqrt(pi)/(2*scl_para))^M)/Gamma(M/2), b=(M/2)-1
a=((gamma(1/2)*(1/(2*sqrt(scl_para))))^nR)/gamma(nR/2);
b=(nR/2)-1;

%Numerical searching of alpha, kappa, eta, mu parameter

%alpha-mu parameter
syms alp mu omeg
A=vpasolve([a==(alp*mu^mu)/(gamma(mu)*omeg^mu), ...
    b==-1+alp*mu],[alp,mu,omeg]);
A

%kappa-mu parameter
% 필요없나?
%gamma(-1)*(1+kap)*mu==exp(-1*kap*mu)*int(fun,t,0,1),
syms kap mu omeg t
fun=exp(kap*mu*t)*(t^mu)*(1-t)^(-2);
fun_int=int(fun,t,0,1);
K=vpasolve([mu==(b+1)/2,
    a==(2*exp(-1*kap*mu)*((1+kap)*mu/omeg)^mu)/gamma(mu)] ...
    ,[kap,mu,omeg]);
K

%eta-mu parameter
syms eta mu omeg
E=vpasolve([mu==(b+1)/4, ...
    (eta+eta^(-1)==((a*gamma(2*mu))/(2*(mu/omeg)^(2*mu)))^(1/mu)-2)] ...
    ,[eta,mu,omeg])