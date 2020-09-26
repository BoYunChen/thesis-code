clear all
close all

dx = 0.0003;
hm = 0.005;
R = 1;
mid = 0.5;
T = 30;

[x,f1,f2,f3] = createf(dx,hm,R,mid);
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
options4 = odeset('Mass',A,'AbsTol',4e-4,'RelTol',3e-2);
options5 = odeset('Mass',A);
[t, y]= ode15s(@PNP2D,[0 T], y,options3);
fprintf('P1 done.\n')

y0 = y(end,:);
[t1, y1] = ode15s(@PNP2D2,[0 0.1],y0,options2);
fprintf('P2 done.\n')
d_his1 = 40*x(0.5*k)*(y1(:,0.5*k+1)-y1(:,0.5*k))*(1+hm)*log(x(0.5*k+1)/x(0.5*k))/dx;


[t2, y2] = ode15s(@PNP2D,[0.1 3.5],y1(end,:),options3);
fprintf('P3 done.\n')
d_his2 = 40*x(0.5*k)*(y2(:,0.5*k+1)-y2(:,0.5*k))*(1+hm)*log(x(0.5*k+1)/x(0.5*k))/dx;

[t3, y3] = ode15s(@PNP2D,[3.5 T],y2(end,:),options3);
fprintf('P4 done.\n')
d_his3 = 40*x(0.5*k)*(y3(:,0.5*k+1)-y3(:,0.5*k))*(1+hm)*log(x(0.5*k+1)/x(0.5*k))/dx;

d_his = [d_his1;d_his2(2:end);d_his3(2:end)];
t_f = [t1;t2(2:end);t3(2:end)];
y_f = [y1;y2(2:end,:);y3(2:end,:)];
y_f(:,1:0.5*k) = y_f(:,1:0.5*k) - d_his;

dlmwrite('x2.csv',x,'precision',10);
dlmwrite('t2.csv',t_f,'precision',10);
dlmwrite('y2.csv',y_f,'precision',10);
