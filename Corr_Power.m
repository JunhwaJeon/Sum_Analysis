clear all; close all;
clc;
%% Correlation Matrix에서 Power Normalizing 되었는지 확인하는 과정
% Uniform correlation Model에서 constraint 이내의 값을 사용하였음
% 잘 되었음을 확인할 수 있음

nT=8; nR=8; %MIMO Scale 지정
sq2 = sqrt(0.5); %상수 지정

phi_T=zeros(nT);
for k=1:nT
    phi_T(k,k)=1;
    for j=k+1:nT
        phi_T(k,j)=-1/(nT-1);
        phi_T(j,k)=-1/(nT-1); %correlation coefficient constraint 범위 내 설정
    end
end
t_corr=sqrtm(phi_T); %Kronecker model correlation mtx.




for i=1:nT
    %norm(t_corr(i,:))^2,norm(t_corr(:,i))^2
    sum=0;
    for j=1:nT
        sum=sum+abs(phi_T(i,j));
    end
    sum^4
end

%{
chann_pow=0;
iter=10^5;
for j=1:iter
    H_w= sq2*(randn(nR,nT)+1j*randn(nR,nT)); %Complex Circular Gaussian channel (Rayleigh)
    H=H_w*t_corr; %channel matrix
    for i=1:nT
        chann_pow=chann_pow+norm(H(i,:))^2;
    end
end
chann_pow/(iter*nT)
%}
