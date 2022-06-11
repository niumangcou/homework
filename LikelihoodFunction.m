clc;clear;
%w:未报告的单位数

%r:报告的故障总数

%t:故障时间
t1=80;%是个区间
t_i1=t1-0.5;
t_i2=t1+0.5;
t_i=(t_i1:0.1:t_i2);
syms r;
syms t;
%B_j:设备寿命
beta=2.788;%beta:形状参数
eta=98;%eta:标度参数
%E_R:威布尔分布的平均值
E_R=eta*gamma(1+1/beta);
%f_T:威布尔分布的概率密度函数
%f_T=beta/eta*(t/eta).^(beta-1)*exp(-(t/eta).^beta);
F=1-exp(-(t/eta).^beta);
f_T=diff(F);
%F_R:退休时间的累积分布公式
F_R=1-exp(-(t/eta).^beta);
%plot(F_R);
FF=beta/eta*(t_i/eta).^(beta-1).*exp(-(t_i/eta).^beta).^2;
%plot(FF);
%P_r:失效时间介于ti1和ti2之间的报废前失效的概率
P_r=int(beta/eta*(r/eta).^(beta-1).*exp(-2*(r/eta).^beta),r,t_i1,t_i2);
%fprintf("P_r=%.5s\n",vpa(rats(P_r)));
zeta=1-int(beta/eta*(r/eta).^(beta-1).*exp(-2*(r/eta).^beta),r,0,118);
%fprintf("zeta=%.5s",vpa(rats(zeta)));

x=diff(f_T*(1-F_R));
a=solve(x,t);
plot(a);
max_A=vpa(max(a));