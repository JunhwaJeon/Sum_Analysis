clear all;
close all;
clc;

sq2 = sqrt(0.5); %상수 지정
nR=4; nT=4;

phi_T=zeros(nT);
p=[0:0.01:1];
A(:,:)=zeros(3,length(p));

for i=1:length(p)
    for k=1:nT
        phi_T(k,k)=1;
        for j=k+1:nT
            phi_T(k,j)=p(i);
            phi_T(j,k)=p(i);
        end
    end
    t_corr=sqrtm(phi_T); %transmit correlation mtx
    sum=0;
    for l=1:nT
        sum=sum+t_corr(l,2);
    end
    scl_para=4*abs(sum)^4;
    a=((gamma(1/2)/(2*sqrt(scl_para)))^nR)/gamma(nR/2);
    b=(nR/2)-1;
    syms alp mu omeg
    [r_alp,r_mu,r_omeg]=vpasolve([a==(alp*mu^mu)/(gamma(mu)*omeg^mu), ...
        b==-1+alp*mu],[alp,mu,omeg],'Random',true);
    %
    if isempty(r_alp)==1 || r_alp==0
        A(1,i)=0;
    else 
        A(1,i)=r_alp;
    end
    %
    if isempty(r_mu)==1 || r_mu==0
        A(2,i)=0;
    else 
        A(2,i)=r_mu;
    end
    %
    if isempty(r_omeg)==1 || r_omeg==0
        A(3,i)=0;
    else 
        A(3,i)=r_omeg;
    end
    %}
end
plot(p,A(1,:),'b-',p,A(2,:),'r-',p,A(3,:),'g-');
xlabel('Uniform Correlation Coefficient');
ylabel('Real part of parameters')
legend('alpha', 'mu', 'omega')
hold on, grid on;