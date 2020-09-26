clear all
close all

dx = 0.0003;
hm = 0.005;
R = 1;
mid = 0.5;
T = 20;

[x,f1,f2,f3] = createf(dx,hm,R,mid); %f1:Na  f2:K  f3:Cl
U = zeros(length(f1),1);

n = 0.3177;
m = 0.05293;
h = 0.5961;

y = [U;f1;f2;f3;m;n;h];
k = (length(y)-3)/4;
A = ones(1,length(y));
A(1:k) = 0*A(1:k);
A = sparse(diag(A));
options1 = odeset('Mass', A,'AbsTol',3e-4,'RelTol',3e-2);
options2 = odeset('Mass',A,'AbsTol',1e-4,'RelTol',1e-2);
options3 = odeset('Mass',A,'AbsTol',2e-4,'RelTol',1e-2);
options4 = odeset('Mass',A,'AbsTol',3e-5,'RelTol',3e-2);

[t, y]= ode15s(@PNP1D, [0 T], y,options1);%3
fprintf('P1 done.\n')

y0 = y(end,:);
[t1, y1] = ode15s(@PNP1D2,[0 0.1],y0,options2);
fprintf('P2 done.\n')
d_his1 = 40*hm*(y1(:,0.5*k+1)-y1(:,0.5*k))/dx;


[t2, y2] = ode15s(@PNP1D,[0.1 3],y1(end,:),options1);
fprintf('P3 done.\n')
d_his2 = 40*hm*(y2(:,0.5*k+1)-y2(:,0.5*k))/dx;


[t3,y3] = ode15s(@PNP1D,[3 70],y2(end,:),options3);
fprintf('P4 done.\n')
d_his3 = 40*hm*(y3(:,0.5*k+1)-y3(:,0.5*k))/dx;

d_his = [d_his1;d_his2(2:end);d_his3(2:end)];
t_f = [t1;t2(2:end);t3(2:end)];
y_f = [y1;y2(2:end,:);y3(2:end,:)];
y_f(:,1:0.5*k) = y_f(:,1:0.5*k) - d_his;

dlmwrite('x1.csv',x,'precision',10);
dlmwrite('t1.csv',t_f,'precision',10);
dlmwrite('y1.csv',y_f,'precision',10);