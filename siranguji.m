clc;%似然估计
clear;
%数据
t_i=[91 41 94 57 32 69 81 94 88 82 32 29 83 87 70 82 106 47 77 66 76 82 71 110 60 51 92 67 107 53 100 76];
r =length(t_i);
A_i=[101 102 102 102 103 104 104 105 106 106 106 106 106 106 106 108 108 111 111 111 112 112 112 113 113 113 114 114 115 116 118 118];
J=14;%批次总数
n_j=[579 1210 598 1223 594 1217 612 1208 1204 616 1208 614 615 589];%安装数量
r_j=[2 1 1 2 3 3 3 2 7 1 2 1 3 1];%报告的故障数量
w_j=n_j-r_j;%未报告的数量
B_j=[118 116 115 114 113 112 111 108 106 105 104 103 102 101];%设备寿命
E_R=98;%报废时间(期望值)
Pr_delta=[0.62 0.31 0.04 0.004 0.004 0.004 0.003 0.003 0.003 0.003 0.001 0.001 0.001 0.001 0.001 0.001];

%参数
[a, b] = wblfit(t_i);%MLE
beta_R=1.5;
eta_R=E_R/gamma(1+1/beta_R);
pi=0;
pi_sum_delta=0;
zeta=0;
zeta_sum_delta=0;
mark_max=-inf;

%公式
%p_r = @(t) (beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta)) .* ( exp(-(t/eta_R).^beta_R) );
p_r2 = @(t) (beta_MLE/eta_MLE.*(t/eta_MLE).^(beta_MLE-1).*exp(-(t/eta_MLE).^beta_MLE)) .* ( exp(-(t/eta_R).^beta_R) );
k_eta=0;

%求pi_i
for eta = b(1,1):10:2000
    k_eta=1+k_eta;
    k_beta=0;
    for beta = b(1,2):0.1: b(2,2)
        k_beta=1+k_beta;
        p_r = @(t) (beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta)) .* ( exp(-(t/eta_R).^beta_R) );
        for i=1:1:r
            t_1=t_i(i)-0.5;
            t_2=t_i(i)+0.5;
            for delta=1:1:length(Pr_delta)
                if t_i(i)+delta-1<A_i(i)
                    pi_i=Pr_delta(delta)*integral(p_r,max(0,t_1),min(t_2,A_i(i)));
                    pi_sum_delta=pi_sum_delta+pi_i;
                end
            end
            pi=log(pi_sum_delta)+pi;
            pi_sum_delta=0;
        end
        for j=1:1:J
            for delta=1:1:length(Pr_delta)
                if t_i(i)+delta-1<B_j(j)
                    
                    zeta_j=Pr_delta(delta)*integral(p_r,0,B_j(j));
                    zeta_sum_delta=zeta_sum_delta+zeta_j;
                end
            end
            zeta_sum_delta=1-zeta_sum_delta;
            zeta=w_j(j)*log(zeta_sum_delta)+zeta;
            zeta_sum_delta=0;
        end
        mark=pi+zeta;
        fprintf("mark=%c\t",mark);

        if mark_max < mark
            fprintf("替换");
            mark_max=max(mark,mark_max);
            eta_MLE=eta;
            beta_MLE=beta;
        end
        pi=0;
        zeta=0;
        fprintf("mark_max=%c\t",mark_max);
        fprintf("eta_max=%c\t",eta_MLE);
        fprintf("beta_max=%c\n",beta_MLE);
    end
end