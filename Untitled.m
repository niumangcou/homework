clc;
clear;
x=(0:1:100);
lamda=0.0092;
E_R=98;
beta_R=1.5;
eta_R=E_R/gamma(1+1/beta_R);
beta=0.5;
eta=590.78;
F=exp(-(x/eta_R).^beta_R);
P=exp(-1/eta_R*x);
p=1/eta_R*exp(-1/eta_R*x);%概率密度
f=diff(F);
h=beta/eta*(x/eta).^(beta-1);
plot(h);
