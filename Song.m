clear;
clc;
data=xlsread('C:\Users\couniumang\Desktop\NewEtaAndBeta.xlsx','Sheet1','A2:C2422');
% data=xlsread('D:NewEtaAndBeta.xlsx','Sheet1','E2:G2422');
eta=data(:,1);
beta=data(:,2);
mark=data(:,3);
m=size(data,1);
xx=unique(data(:,1));
yy=unique(data(:,2));
[x,y]=meshgrid(xx,yy);
z=0*x.*y;
for i=1:m
    k1=find(xx==data(i,1));
    k2=find(yy==data(i,2));
    z(k2,k1)=z(k2,k1)+1;
end
%  bar3(z)
%  set(gca,'xticklabel',xx,'yticklabel',yy)
k=1;
for j=1:1:length(data)
    if (beta(j)>=3 && beta(j)<=4 && eta(j)>=300 && eta(j)<=400)%163 1838 (33.6,131.9)
        BBB(k,1)=eta(j);
        BBB(k,2)=beta(j);
        BBB(k,3)=mark(j);
        k=k+1;
    end
end