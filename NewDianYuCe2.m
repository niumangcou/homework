%区间预测
clc;clear;
t_i = [40 9 28 17 23 36 13 23 14 33 15 36 22 46 7 18 6 23 24 27 23 46 34 38 21 30 14 29 15 16 39 12 7 34 14 22 31 34 34 26 32 26 9 34 45 35 8 20 18 26 22 15 33 27 17 48 34 49 23 49 19 28 21 34 14 45 17 22 15 28];
r = length(t_i);
n_j = 10000;
w_j=n_j-r;%未报告的数量
rho_j(1)=0.0086;
n_xing=w_j;
w=2*pi/(n_xing+1);

%由F_N (n),n=0,1,…,n^*表示的N(s)的cdf
sum_B=0;
A=1/(n_xing+1);
for n_P=1:1:n_xing
    fprintf("n_P=%d\n",n_P);
for l=1:1:n_xing
    prod_C=1;
     %fprintf("l=%d\n",l);
    B=(exp(-1i * w *l * n_P)-exp(-1i * w * l))/(1-exp(-1i * w * l));
       C=(1 - rho_j(1) + rho_j(1) * exp(1i * w * l)) ^ w_j;
       prod_C=prod_C * C;
       %fprintf("prod_C:%c\n",prod_C);
   sum_B=sum_B + B * prod_C;
   %fprintf("sum_B=%s\n",sum_B);
end
FN_n(n_P)=A * sum_B;
end
plot(FN_n);