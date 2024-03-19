function moment=SumMoment_4(al, om, n)
% al (shape parameter) 
% om (scale parameter) 로 갖는Line 3 % om (scale parameter) 로 갖는 summand
% 4개를 더했을 때의 sum의 n차 모멘트

N=[1:1:n];
summand_moment=zeros(length(N)+1);

% Summand의 n차 모멘트 구하기
for i=1:length(N)+1
    if i==1
        summand_moment(i)=1;
    else
    summand_moment(i)=om^(N(i-1)/al)*gamma(1+N(i-1)/al);
    end
end
moment_tmp=1;

% Sum moment 루프
for n1=0:n
    for n2=0:n1
        for n3=0:n2
            moment_tmp=moment_tmp*nchoosek(n2,n3).*summand_moment(n3+1);
        end
        moment_tmp=moment_tmp*nchoosek(n1,n2)*summand_moment(n2+1);
    end
    moment_tmp=moment_tmp*nchoosek(n,n1)*summand_moment(n1+1);
end
moment=moment_tmp;
end