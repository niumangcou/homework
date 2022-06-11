clc;clear;
t_i = [40 30 28 17 23 36 42 23 40 33 15 36 47 46 31 18 6 23 24 27 23 46 34 38 32 30 27 47 33 16 39 26 12 34 31 22 31 32 38 41 32 41 29 38 45 35 8 20 49 41 45 46 33 42 40 48 38 40 39 49 19 28 21 34 25 45 17 22 33 28];
r = length(t_i);
A_i = 50;
n_j = 10000;
w_j = n_j - r;
E_R = 36;
Pr_delta=[0.68 0.24 0.04 0.01 0.01 0.005 0.005 0.005 0.001 0.001 0.001 0.001 0.001];
%参数
[a, b] = wblfit(t_i);
beta_R = 1.5;
eta_R = E_R / gamma( 1 + 1 / beta_R );
pi = 0;
pi_sum_delta = 0;
zeta = 0;
zeta_sum_delta = 0;
N=200;
for n=1:1:N
    fprintf("n=%d\n",n);
    mark_max = - inf;
    for eta = b(1,1):10:1000
        for beta = b(1,2)*0.3:0.1: 5
            p_r = @(t) (beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta)) .* ( exp(-(t/eta_R).^beta_R) );
            for i=1:1:r
                W_i=frnd(10,50,r,1);%随机权重
                t_1=t_i(i)-0.5;
                t_2=t_i(i)+0.5;
                pi_i=integral(p_r,max(0,t_1),min(t_2,A_i));
                for delta=1:1:length(Pr_delta)
                    
                    if t_i(i)+delta-1<A_i
                        pi_delta=Pr_delta(delta)*pi_i;
                        pi_sum_delta=pi_sum_delta+pi_delta;
                    end
                end
                pi=W_i(i,1)*log(pi_sum_delta)+pi;%*
                pi_sum_delta=0;
            end
            W_jk=frnd(10,50,w_j,1);%随机权重
            zeta_j=integral(p_r,0,A_i);
            for delta=1:1:length(Pr_delta)
                if t_i(i)+delta-1<A_i
                    zeta_delta=Pr_delta(delta)*zeta_j;
                    zeta_sum_delta=zeta_sum_delta+zeta_delta;
                end
            end
            zeta_sum_delta=1-zeta_sum_delta;            
            zeta=sum(W_jk)*log(zeta_sum_delta)+zeta;%*
            zeta_sum_delta=0;
            mark=pi+zeta;
            if mark_max < mark
                mark_max=max(mark,mark_max);
                AAA(n,3)=mark_max;
                AAA(n,1)=eta;
                AAA(n,2)=beta;
            end
            pi=0;
            zeta=0;
        end
    end
end
% fitdist(nn,'Normal');