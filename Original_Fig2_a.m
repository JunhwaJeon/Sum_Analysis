clear all, close all;
clc;
%% 원 논문 Fig2_a
%% Transmit Antenna Selection Gain 영향 분석
%% Quantization bit=3, Nr=1, Nt=1,4,16,64 에서의 Transmit SNR - Ergodic Rate 그래프
%% 파라미터 설정
T_SNR_dB=[-10:1:20]; %SNR 범위 설정
T_SNR_linear=10.^(T_SNR_dB/10); %linear 스케일 SNR설정
N_iter=5000; %반복 횟수 (Ergodic capacity 구하기 위해서) 
sq2 = sqrt(0.5); %상수 지정
nT=0; nR=1; %MIMO Scale 지정
q_gain=0.96546;

for Icase=1:4 %Tx개수에 따른 그래프
    if Icase==1, nT=1; 
    elseif Icase==2, nT=4; 
    elseif Icase==3, nT=16; 
    else nT=64; 
    end

    R(Icase,:) = zeros(1,length(T_SNR_dB));%Capacity 정보 담을 행렬 지정 (안테나 경우*SNR 범위)
    R_candi=linspace(0,0,nT); %Maximum 선택 위한 후보값 담을 행렬(벡터) 지정

    %% Ergodic Rate 계산
    for i=1:length(T_SNR_dB)
        for iter=1:N_iter %반복
            H= sq2*(randn(nR,nT)+1j*randn(nR,nT)); %Complex Circular Gaussian channel (Rayleigh)
            for j=1:nT
                sum_four_sqr=0;
                norm_sqr=norm(H(:,j))^2;
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

plot(T_SNR_dB,R(1,:),'k--', T_SNR_dB,R(2,:),'b-', T_SNR_dB,R(3,:),'b-');
hold on, grid on,
plot(T_SNR_dB,R(4,:),'b-');
xlabel('Transmit SNR[dB]'); ylabel('Ergodic Rate [bps/Hz]');