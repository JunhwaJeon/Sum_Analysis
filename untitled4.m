% 파라미터 설정
lambda = 64;
k = 1/2;

% 범위 설정
x = 0:0.1:200;

% Weibull 분포의 확률 밀도 함수(PDF) 계산
pdf = (k/lambda) * (x/lambda).^(k-1) .* exp(-(x/lambda).^k);

% 그래프 그리기
plot(x, pdf, 'LineWidth', 2);
title('Weibull Distribution (scale=64, shape=1/2)');
xlabel('x');
ylabel('PDF');
grid on;