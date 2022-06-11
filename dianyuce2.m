%区间预测
clc;clear;
n_j=[579 1210 598 1223 594 1217 612 1208 1204 616 1208 614 615 589];%安装数量
r_j=[2 1 1 2 3 3 3 2 7 1 2 1 3 1];%报告的故障数量
w_j=n_j-r_j;%未报告的数量
rho_j=[4.03707918689982e-08 1.17629241739095e-07 2.64681441072803e-07 5.04457488693626e-07 8.59497817971477e-07 1.35195171813679e-06 2.00356957304952e-06 2.83569269816493e-06 3.86924238202847e-06 5.12470888741655e-06 6.62214079367732e-06 8.38113487990766e-06 1.04208266533107e-05 1.27598815744633e-05 1.54164870010438e-05 1.84083448537943e-05 2.17526649978593e-05 2.54661593262441e-05 2.95650365283233e-05 3.40649975241190e-05 3.89812315438656e-05 4.43284128318151e-05 5.01206979530867e-05 5.63717236824691e-05 6.30946054543546e-05 7.03019363533552e-05 7.80057866255792e-05 8.62177036910091e-05 9.49487126378959e-05 0.000104209317180564 0.000114009501062491 0.000124358729886997 0.000135265953358321 0.000146739607916304 0.000158787619748333 0.000171417408162591 0.000184635889307073 0.000198449480219209 0.000212864103191294 0.000227885190437312 0.000243517689047062 0.000259766066213860 0.000276634314722404 0.000294125958683713 0.000312244059504391 0.000330991222077721 0.000350369601184462 0.000370380908091447 0.000391026417336429 0.000412306973687855 0.000434222999268556 0.000456774500832602 0.000479961077184825 0.000503781926732814 0.000528235855161392 0.000553321283219877 0.000579036254612674 0.000605378443983976 0.000632345164987605 0.000659933378433265 0.000688139700500714 0.000716960411013588 0.000746391461764843 0.000776428484886012 0.000807066801252695 0.000838301428918900 0.000870127091573115 0.000902538227009146 0.000935528995605020 0.000969093288803429 0.00100322473758739 0.00103791672094505 0.00107316237431764 0.00110895459802504 0.00114528606566312 0.00118214923246796 0.00121953634364133 0.00125743944263290 0.00129585037937405 0.00133476081845894 0.00137416224726826 0.00141404598403150 0.00145440318582355 0.00149522485649190 0.00153650185451050 0.00157822490075676 0.00162038196634592 0.00166296585863527 0.00170596691616933 0.00174937536986378 0.00179318135036382 0.00183737489530781 0.00188193887286981 0.00192686954995580 0.00197215669453763 0.00201779000917507 0.00206375487286990 0.00211004475274863 0.00215664917818719 0.00220338244542899];
n_xing=sum(w_j);
w=2*pi/(n_xing+1);
J=14;
%由F_N (n),n=0,1,…,n^*表示的N(s)的cdf
sum_B=0;
A=1/(n_xing+1);
for n_P=1:1:n_xing
for l=1:1:n_xing
    prod_C=1;
     fprintf("l=%d\n",l);
    B=(exp(-1i * w *l * n_P)-exp(-1i * w * l))/(1-exp(-1i * w * l));
   for j=1:1:J
       C=(1 - rho_j(1) + rho_j(1) * exp(1i * w * l)) ^ w_j(j);
       prod_C=prod_C * C;
       %fprintf("prod_C:%c\n",prod_C);
   end
   sum_B=sum_B + B * prod_C;
   fprintf("sum_B=%s\n",sum_B);
end
FN_n(n_P)=A * sum_B;
end
plot(FN_n);