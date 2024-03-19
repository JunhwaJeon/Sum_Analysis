clear all
close all
clc
% 파라미터 설정
lambda = 64;
k = 1/2;
num_variables = 4; % 합할 확률 변수의 개수

% 범위 설정
x = 0:0.1:1000;

% 각 확률 변수의 PDF 계산
pdf_each = (k/lambda) * (x).^(k-1) .* exp((-(x).^k)/lambda);

% 확률 변수들의 합의 PDF 계산
pdf_sum = pdf_each;
for i = 2:num_variables
    pdf_sum = conv(pdf_sum, pdf_each);
end

% Moment Matching 그래프
[al, mu, om]=Weibull_MoM_4(k, lambda);
[al, mu, om]
app_al_mu=(al*(mu^mu)*(x.^(al*mu-1))).*exp((-1*mu*(x.^al))/om)/(gamma(mu)*(om^mu));

% 그래프 그리기
plot(pdf_sum(30001:40001), 'LineWidth', 2);
hold on;
plot(app_al_mu(1:10001), 'r-');
title('PDF of Sum of 4 Weibull Distributed Variables');
xlabel('x');
ylabel('PDF');
grid on;