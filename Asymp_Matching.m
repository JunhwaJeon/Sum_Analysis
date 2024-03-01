clear all;
close all;
clc;

%% Channel Matrix Under Spatial Correlation Fromulation%%
sq2 = sqrt(0.5); %상수 지정

%iid Rayleigh channel matrix
nR=8,nT=8;
H_w=sq2*(randn(nR,nT)+1j*randn(nR,nT));

%spatial correlation matrix
%random correlation coefficient 주어서 어떤 distribution이 matching이 잘 되는지 확인
phi_T=zeros(nT);
for i=1:nT
    phi_T(i,i)=1;
    for j=i+1:nT
        x=-1+2*rand;
        phi_T(i,j)=x+1j*sqrt(1-x^2);
        phi_T(j,i)=conj(phi_T(i,j));
    end
end
t_corr=sqrtm(phi_T);

%channel matrix
H=H_w*t_corr;

%% Define Exact Sum %%

%scale parameter 구하기
sum=0;
for i=1:nT
    sum=sum+H(i,2);
end
scl_para=sum^4;

%Weibull distribution 정의하기
%pdf_wei(x)=(1/(scl_para*2))*((x/scl_para)^(-1/2))*exp(-(x/scl_para)^(1/2));
%Weibull sum 정의하기 - 다중적분 이용하여 정의
%sum_wei(x)=

%% Define Original Approximation %%

%% Define Asymptotic Matching Approximation %%

%Summand's coefficient are given by: a=1/(2k), b=-1/2
%Sum coefficient are given by: a=((sqrt(pi)/(2*scl_para))^M)/Gamma(M/2), b=(M/2)-1
a=((sqrt(pi)/(2*scl_para))^nR)/gamma(nR/2);
b=(nR/2)-1;

%Numerical searching of alpha, kappa, eta, mu parameter
%alpha-mu distribution
syms alp mu omeg
a=(alp*mu^mu)/(gamma(mu)*omeg^mu);
X=vpaslove([a==(alp*mu^mu)/(gamma(mu)*omeg^mu) ...
    b==-1+alp*mu],[alp,mu,omeg]);
X