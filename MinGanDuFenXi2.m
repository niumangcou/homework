clc;clear;%点预测

eta_MLE=[[193.286761946048 233.286761946048 253.286761946048]];
beta_R=[1 1.5 2];
beta_MLE=[[2.89504932055273 2.69504932055273 2.59504932055273]];
%数据
t_i = [40 30 28 17 23 36 42 23 40 33 15 36 47 46 31 18 6 23 24 27 23 46 34 38 32 30 27 47 33 16 39 26 12 34 31 22 31 32 38 41 32 41 29 38 45 35 8 20 49 41 45 46 33 42 40 48 38 40 39 49 19 28 21 34 25 45 17 22 33 28];
r = length(t_i);
A_i = 50;
n_j = 10000;
w_j = n_j - r;
E_R = 36;
Pr_delta=[0.68 0.24 0.04 0.01 0.01 0.005 0.005 0.005 0.001 0.001 0.001 0.001 0.001];


s=1;%几个月内
for k=1:1:length(beta_R)
    eta_R = E_R / gamma( 1 + 1 / beta_R(k) );
    p_r2 = @(t) (beta_MLE(k)/eta_MLE(k).*(t/eta_MLE(k)).^(beta_MLE(k)-1).*exp(-(t/eta_MLE(k)).^beta_MLE(k))) .* ( exp(-(t/eta_R).^beta_R(k)) );
for point=1:1:100

    rho_j(point,k)=0;
    for l=0:1:s-1
    gamma_sum_delta=0;
    xi_sum_delta=0;
%         if point<B_j(j)% 原来是R
            B_jl1=max(point,point+l-0.5);
            B_jl2=point+l+0.5;
            for delta=1:1:length(Pr_delta) %加入延迟报告
                if (B_jl1<=point+delta-1 && point+delta-1<=B_jl2)
                    fprintf("%d\t%d\t%c\n",point+delta-1,B_jl1,B_jl2);
                    gamma_j=Pr_delta(delta)*integral(p_r2,0,B_jl2);
                    gamma_sum_delta=gamma_sum_delta+gamma_j;
                end
                if (0<=point+delta-1 && point+delta-1<=A_i)
                    xi_j=Pr_delta(delta)*integral(p_r2,0,point);
                    xi_sum_delta=xi_sum_delta+xi_j;
                end
            end
            xi_sum_delta=1-xi_sum_delta;
            h_j=gamma_sum_delta/xi_sum_delta;
        rho_j(point,k)=rho_j(point,k)+h_j; %批次j中的单元在DFD后s个月内被报告为故障的概率
        %gamma_sum_delta=0;
        fprintf("%d\n",xi_sum_delta);
        
%         end
    end
    n(point,k)=rho_j(point,k)*w_j;%点预测
end

end
y=plot(n),legend('E_R=30', 'E_R=36', 'E_R=40');%分布曲线
xlabel('时间点');
ylabel('故障数');