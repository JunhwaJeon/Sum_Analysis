close all; clear; clc;

%% weibull 분포 k=1/2 omega=1 일 때 moment 
% n차 모멘트는 gamma(1+2n)
W_moment=[gamma(3),gamma(5),gamma(7),gamma(9),gamma(11)];


%
function result=sum_moment(n,n_R)
for n1=0:1:n
    r=r*nchoosek(n,n1)*gamma(1+2*(n-n1));
end
end

function total_sum=recursive_sum(n, n_R, W_expectations, )
end