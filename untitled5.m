clear all; close; clc;

al=1/2; om=64;
fir_moment=SumMoment_4(al, om, 1);
sec_moment=SumMoment_4(al, om, 2);
four_moment=SumMoment_4(al, om, 4);
z1=log((fir_moment)^2/(sec_moment));
z2=log((sec_moment)^2/(four_moment));

% 이변수 함수 정의
syms m k
fun=equa(m,k,z1,z2);

% 초기 추정값 설정
x0 = [1.0 1.0];

% 3. 최적화 함수 호출
options = optimoptions('fminunc', 'Display', 'iter'); % 옵션 설정
[x_min, fval, exitflag, output] = fminunc(fun, x0, options); % 최적화 함수 호출

% 4. 최적화 결과 확인
disp(['최소값: ', num2str(fval)]);
disp(['최적해: ', num2str(x_min)]);
disp(['종료 상태: ', num2str(exitflag)]);
disp(['출력 정보:']);
disp(output);

% 결과 출력
disp('Optimal point:');
disp(x);
disp('Function value at optimal point:');
disp(fval);

function F=equa(m,k,z1,z2)
F=(g(m,k)-log(z1))^2+(g(m,k/2)-log(z2))^2;
end

function result = g(m, k)
    h1 = 2*gammaln(m + 1/k);
    h2 = gammaln(m);
    h3 = gammaln(m + 2/k);
    
    result = h1 - h2 - h3;
end