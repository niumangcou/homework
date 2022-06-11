clc;clear;%点预测
eta_R=108.5578;
%eta_MLE=89.0305;
% eta_MLE=1394;
eta_MLE=644;
beta_R=1.5;
%beta_MLE=3.3876;
% beta_MLE=2.98;
beta_MLE=2.98;
%数据
t_i=[91 41 94 57 32 69 81 94 88 82 32 29 83 87 70 82 106 47 77 66 76 82 71 110 60 51 92 67 107 53 100 76];%32
r =length(t_i);
A_i=[101 102 102 102 103 104 104 105 106 106 106 106 106 106 106 108 108 111 111 111 112 112 112 113 113 113 114 114 115 116 118 118];%32
J=14;%批次总数
%n_j=[5795 12100 5985 12233 5946 12175 6124 12083 12040 6166 12080 6147 6155 5892];%安装数量
n_j=[579 1210 598 1223 594 1217 612 1208 1204 616 1208 614 615 589];%安装数量
r_j=[2 1 1 2 3 3 3 2 7 1 2 1 3 1];%报告的故障数量
w_j=n_j-r_j;%未报告的数量
B_j=[118 116 115 114 113 112 111 108 106 105 104 103 102 101];%设备寿命14
R=98;%退休时间
E_R=98;%报废时间(期望值)
Pr_delta=[0.62 0.31 0.04 0.004 0.004 0.004 0.003 0.003 0.003 0.003 0.001 0.001 0.001 0.001 0.001 0.001];

p_r2 = @(t) (beta_MLE/eta_MLE.*(t/eta_MLE).^(beta_MLE-1).*exp(-(t/eta_MLE).^beta_MLE)) .* ( exp(-(t/eta_R).^beta_R) );

s=5;%几个月内


for point=1:1:300
for j=1:1:J
    rho_j(j,point)=0;
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
                if (0<=point+delta-1 && point+delta-1<=B_j(j))
                    xi_j=Pr_delta(delta)*integral(p_r2,0,point);
                    xi_sum_delta=xi_sum_delta+xi_j;
                end
            end
            xi_sum_delta=1-xi_sum_delta;
            h_j=gamma_sum_delta/xi_sum_delta;
        rho_j(j,point)=rho_j(j,point)+h_j; %批次j中的单元在DFD后s个月内被报告为故障的概率
        %gamma_sum_delta=0;
        fprintf("%d\n",xi_sum_delta);
        
%         end
    end
    n(j,point)=rho_j(j,point)*w_j(j);%点预测
end
end
y=plot(n(1,:));%分布曲线