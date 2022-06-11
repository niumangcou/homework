clc;clear;
t_i=[91 41 94 57 32 69 81 94 88 82 32 29 83 87 70 82 106 47 77 66 76 82 71 110 60 51 92 67 107 53 100 76];
y=plot(ksdensity(t_i));%分布曲线
m=ksdensity(t_i,0.1,'function','pdf');
t_0001=quantile(t_i,0.001);