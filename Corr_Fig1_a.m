clear all, close all
clc;

%{
Correlation Matrix는 Uniform Correlation Model
%}

%% 파라미터 설정
T_SNR_dB=[-15:1:20]; %SNR 범위 설정
T_SNR_linear=10.^(T_SNR_dB/10); %linear 스케일 SNR설정
N_iter=1000; %반복 횟수 (Ergodic capacity 구하기 위해서) 
sq2 = sqrt(0.5); %상수 지정
nT=32; nR=8; %MIMO Scale 지정
q_gain=0;

phi_T=zeros(nT);
for k=1:nT
    phi_T(k,k)=1;
    for j=k+1:nT
        phi_T(k,j)=-1/(nT-1);
        phi_T(j,k)=-1/(nT-1); %correlation coefficient=1/2
    end
end
t_corr=sqrtm(phi_T); %Kronecker model correlation mtx.


%% Quantization bit 지정, b에 따른 상수 지정, b infty인 경우 근사식 이용
for Icase=1:5 %quantization bit에 따른 그래프
    if Icase==1, q_gain=0.6364; 
    elseif Icase==2, q_gain=0.8825; 
    elseif Icase==3, q_gain=0.96546; 
    elseif Icase==4, q_gain=0.990503; 
    else q_gain=1; 
    end

n=min(nT,nR); %channel mtx 짧은 부분 크기
I = eye(n); %해당 스케일의 단위행렬 지정

R(Icase,:) = zeros(1,length(T_SNR_dB));%Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)
R_candi=linspace(0,0,nT); %Maximum 선택 위한 후보값 담을 행렬(벡터) 지정

    %% Ergodic Capacity 계산
    for i=1:length(T_SNR_dB)
        for iter=1:N_iter %반복
            H_w= sq2*(randn(nR,nT)+1j*randn(nR,nT)); %Complex Circular Gaussian channel (Rayleigh)
            H=H_w*t_corr; %channel matrix
         if nR>=nT, HH = H'*H; else HH = H*H'; end
            for j=1:nT
            sum_four_sqr=0; norm_sqr=0;
            norm_sqr=(H(:,j))'*H(:,j);
                for k=1:nR
                sum_four_sqr=sum_four_sqr+abs(H(k,j))^4;
                end
                R_candi(j)=log2(1+(T_SNR_linear(i).*q_gain.*(norm_sqr).^2)/(norm_sqr+T_SNR_linear(i).*(1-q_gain).*sum_four_sqr));
            end
            R(Icase,i)=R(Icase,i)+max(R_candi);
        end
    end
end

R = R/N_iter; %Expectation 계산

plot(T_SNR_dB,R(1,:),'b-', T_SNR_dB,R(2,:),'b-', T_SNR_dB,R(3,:),'b-');
hold on, grid on,
plot(T_SNR_dB,R(4,:),'b-', T_SNR_dB, R(5,:),'b-');
xlabel('Transmit SNR[dB]'); ylabel('Ergodic Rate [bps/Hz]');