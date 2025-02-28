clear all, close all;
clc;
%% 원 논문 Fig2_b
%% Transmit Antenna Selection Gain 영향 분석
%% Transmit SNR=5dB, Quantization bit=2,3,4 에서의 Tx-Ergodic Rate 그래프
%% 파라미터 설정
T_SNR_dB=5; %Transmit SNR
T_SNR_linear=10.^(T_SNR_dB/10); %linear 스케일 SNR설정
N_iter=10^5; %반복 횟수 (Ergodic capacity 구하기 위해서) 
sq2 = sqrt(0.5); %상수 지정
nT=[1:3:64]; 
nR=1; %MIMO Scale 지정

%% Quantization bit 지정, b에 따른 상수 지정, b infty인 경우 근사식 이용
for Icase=1:3 %quantization bit에 따른 그래프
    if Icase==1, q_gain=0.8825; 
    elseif Icase==2, q_gain=0.96546; 
    else q_gain=0.990503; 
    end

    R(Icase,:) = zeros(1,length(nT));%Capacity 정보 담을 행렬 지정 (Qbit 수*nT범위)

    %% Ergodic Capacity 계산
    for i=1:length(nT)
        R_candi=linspace(0,0,nT(i)); %Maximum 선택 위한 후보값 담을 행렬(벡터) 지정
        for iter=1:N_iter %반복
            H= sq2*(randn(nR,nT(i))+1j*randn(nR,nT(i))); %Complex Circular Gaussian channel (Rayleigh)
            for j=1:nT(i)
                sum_four_sqr=0;
                norm_sqr=norm(H(:,j))^2;
                for k=1:nR
                    sum_four_sqr=sum_four_sqr+abs(H(k,j))^4;
                end
                R_candi(j)=log2(1+(T_SNR_linear*q_gain*(norm_sqr)^2)/(norm_sqr+T_SNR_linear*(1-q_gain)*sum_four_sqr));
            end
            R(Icase,i)=R(Icase,i)+max(R_candi);
        end
    end
end
R = R/N_iter; %Expectation 계산



plot(nT,R(1,:),'b-', nT,R(2,:),'b-', nT,R(3,:),'b-');
hold on, grid on,
xlim([1,64]);
xlabel('Number of Transmit Antennas Nt'); ylabel('Ergodic Rate [bps/Hz]');