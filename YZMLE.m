clc;%似然估计
clear;
%数据
t_i=[16 39 50 22 51 31 23 46 41 51 50 44 52 49 53 52 42 51 49 53 24 63 28 43 23 57 43 43 60 56 61 41 57 53 66 56 47 34 53 30 60 28 67 43 60 54 43 33];
r =length(t_i);
A_i=[57 57 58 58 58 59 59 59 59 60 60 60 61 61 62 62 62 62 63 63 64 64 64 65 65 65 66 66 66 66 67 67 67 67 67 67 68 68 68 68 69 69 69 70 70 70 70 70];
J=14;%批次总数
n_j=[1579 2029 2610 4722 2440 829 830 1669 3348 3475 5128 5003 1048 840];%安装数量
r_j=[2 3 4 3 2 4 2 3 3 4 6 4 3 5];%报告的故障数量
w_j=n_j-r_j;%未报告的数量
B_j=[70 69 68 67 66 65 64 63 62 61 60 59 58 57];%设备寿命
E_R=50;%报废时间(期望值)
Pr_delta=[0.68 0.24 0.04 0.01 0.01 0.005 0.005 0.005 0.001 0.001 0.001 0.001 0.001];

%参数
[a, b] = wblfit(t_i);%MLE
beta_R=1.5;
eta_R=E_R/gamma(1+1/beta_R);
pi=0;
zeta=0;
zeta_sum_delta=0;
mark_max=-inf;

for beta = b(1,2)*0.1:0.1: 4
    for eta = b(1,1):10:2000
        p_r = @(t) (beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta)) .* ( exp(-(t/eta_R).^beta_R) );
        for i=1:1:r
            t_1=t_i(i)-0.5;
            t_2=t_i(i)+0.5;
            pi_i=integral(p_r,max(0,t_1),min(t_2,A_i(i)));
            pi_sum_delta=0;
            for delta=1:1:length(Pr_delta)
                if t_i(i)+delta-1<A_i(i)
                    pi_delta=Pr_delta(delta)*pi_i;
                    pi_sum_delta=pi_sum_delta+pi_delta;
                end
            end
            pi=log(pi_sum_delta)+pi;
            
        end
        
        for j=1:1:J
            zeta_j=integral(p_r,0,B_j(j));
            for delta=1:1:length(Pr_delta)
                if t_i(i)+delta-1<B_j(j)
                    zeta_delta=Pr_delta(delta)*zeta_j;
                    zeta_sum_delta=zeta_sum_delta+zeta_delta;
                end
            end
            zeta_sum_delta=1-zeta_sum_delta;
            zeta=w_j(j)*log(zeta_sum_delta)+zeta;
            zeta_sum_delta=0;
        end
        mark=pi+zeta;
        if mark_max < mark
            mark_max=max(mark,mark_max);
            eta_MLE=eta;
            beta_MLE=beta;
            fprintf("mark_max=%c\t",mark_max);
            fprintf("eta_max=%c\t",eta_MLE);
            fprintf("beta_max=%c\n",beta_MLE);
        end
        pi=0;
        zeta=0;
    end
end%356.95 3.35