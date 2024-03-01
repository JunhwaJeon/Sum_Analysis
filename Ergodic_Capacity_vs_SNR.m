% Ergodic_Capacity_vs_SNR.m
clear all, close all

%%파라미터 설정
SNR_dB=[0:5:20]; %SNR 범위 설정
SNR_linear=10.^(SNR_dB/10); %그래프 SNR 스케일 설정
N_iter=1000; %반복 횟수 (Ergodic capacity 구하기 위해서) 
sq2 = sqrt(0.5); %상수 지정

%%MIMO 안테나 개수 지정
for Icase=1:5
if Icase==1, nT=1; nR=1; % 1x1
    elseif Icase==2, nT=1; nR=2; % 1x2
elseif Icase==3, nT=2; nR=1; % 2x1
elseif Icase==4, nT=2; nR=2; % 2x2
else nT=4; nR=4; % 4x4
end

n=min(nT,nR); %channel mtx 짧은 부분 크기
I = eye(n); %해당 스케일의 단위행렬 지정

C(Icase,:) = zeros(1,length(SNR_dB));%Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)

%% Ergodic Capacity 계산
for iter=1:N_iter %반복
H = sq2*(randn(nR,nT)+j*randn(nR,nT)); %Complex Circular Gaussian channel (Rayleigh)
if nR>=nT, HH = H'*H; else HH = H*H'; end
for i=1:length(SNR_dB) % Random channel generation
C(Icase,i) = C(Icase,i)+log2(real(det(I+SNR_linear(i)/nT*HH)));
end
end
end
C = C/N_iter;
plot(SNR_dB,C(1,:),’b-o’, SNR_dB,C(2,:),’b-’, SNR_dB,C(3,:),’b-s’);
hold on, plot(SNR_dB,C(4,:),’b->’, SNR_dB,C(5,:),’b-^’);
xlabel(’SNR[dB]’); ylabel(’bps/Hz’);