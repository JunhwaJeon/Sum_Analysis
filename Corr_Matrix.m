%% Correlation Matrix의 생성 확인
%% Uniform Correlaiton Model에서 Nt 달라짐에 따라 바뀌는 Correlation Matrix 형태 확인

nT=[4 8 16]; nR=8; %MIMO Scale 지정
sq2 = sqrt(0.5); %상수 지정
A(:)=linspace(0,0,length(nT));

for scale=1:length(nT)
    phi_T=zeros(nT(scale));
    for k=1:nT(scale)
        phi_T(k,k)=1;
        for j=k+1:nT(scale)
            %phi_T(k,j)=-1/(nT(scale)-1);
            %phi_T(j,k)=-1/(nT(scale)-1); %correlation coefficient constraint 범위 내 설정
            phi_T(k,j)=1/(2);
            phi_T(j,k)=1/(2);
        end
    end
    t_corr=sqrtm(phi_T); %Kronecker model correlation mtx.
    t_corr
end
