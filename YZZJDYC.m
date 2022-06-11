%clc;clear;%点预测
beta_R=1.5;
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

eta_R=E_R/gamma(1+1/beta_R);
%data=xlsread('C:\Users\couniumang\Desktop\NewEtaAndBeta.xlsx','Sheet1','A2:C1198');
s=1;%几个月内
for i=1:1:length(BBB)
for j=1:1:J
    eta_MLE=BBB(i,1);
    beta_MLE=BBB(i,2);
    p_r2 = @(t) (beta_MLE/eta_MLE.*(t/eta_MLE).^(beta_MLE-1).*exp(-(t/eta_MLE).^beta_MLE)) .* ( exp(-(t/eta_R).^beta_R) );
for point=1:1:150
    
    rho_j(i,point)=0;
    for l=0:1:s-1
    gamma_sum_delta=0;
    xi_sum_delta=0;
%         if point<B_j(j)% 原来是R
            B_jl1=max(point,point+l-0.5);
            B_jl2=point+l+0.5;
            gamma_j=integral(p_r2,0,B_jl2);
            xi_j=integral(p_r2,0,point);
            for delta=1:1:length(Pr_delta) %加入延迟报告
                if (B_jl1<=point+delta-1 && point+delta-1<=B_jl2)
                    fprintf("%d\t%d\t%c\n",point+delta-1,B_jl1,B_jl2);
                    gamma_delta=Pr_delta(delta)*gamma_j;
                    gamma_sum_delta=gamma_sum_delta+gamma_delta;
                end
                if (0<=point+delta-1 && point+delta-1<=B_j(j))
                    xi_delta=Pr_delta(delta)*xi_j;
                    xi_sum_delta=xi_sum_delta+xi_delta;
                end
            end
            xi_sum_delta=1-xi_sum_delta;
            h_j=gamma_sum_delta/xi_sum_delta;
        rho_j(i,point)=rho_j(i,point)+h_j; %批次j中的单元在DFD后s个月内被报告为故障的概率
        fprintf("%d\n",xi_sum_delta);
%         end
    end
    n(i,point)=rho_j(i,point)*sum(w_j);%点预测
end
end
end
for i=1:1:length(BBB)
   y=plot(n(i,:));
   hold on
end
% y=plot(n),legend('s=1', 's=2', 's=3');%分布曲线
% xlabel('时间点');
% ylabel('故障数');