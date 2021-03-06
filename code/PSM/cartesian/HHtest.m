clear all
close all
format long

T = 25;

m = 0.0566;
n = 0.3248;
h = 0.5778;
Vm = -71.4;

y = [m;n;h;Vm];
[t_HH1, y_HH1]= ode15s(@HH2, [0 0.1], y);
[t_HH2, y_HH2]= ode15s(@HH1, [0.1 T], y_HH1(end,:));

t_HH = [t_HH1(1:end-1);t_HH2];
y_HH = [y_HH1(1:end-1,:);y_HH2];
dlmwrite('t_HH.csv',t_HH,'precision',10)
dlmwrite('y_HH.csv',y_HH,'precision',10)
%plot(t_HH,y_HH(:,end))
%--------------------------------------------------------
% figure(1)
% plot(t_HH,y_HH(:,end),'-*')
% grid on
% figure(2)
% plot(t_HH,y_HH(:,end-3:end-1),'-*')
% legend('m','n','h')
% grid on
%---------------------------------------------------------
% dx = 0.0004;
% hm = 0.005;
% R = 1;
% mid = 0.5;
% 
% [x,f1,f2,f3] = createf(dx,hm,R,mid);
% U = zeros(length(f1),1);
% 
% y = [U;f1;f2;f3;m;n;h];
% k = (length(y)-3)/4;
% A = ones(1,length(y));
% A(1:k) = 0*A(1:k);
% A = sparse(diag(A));
% options = odeset('Mass', A,'AbsTol',0.0002);
% [t_PNP, y]= ode15s(@PNP1D, [0 T], y,options);
% d_his = -24.1*40*hm*(y(:,0.5*k+1)-y(:,0.5*k))/dx;
% 
% %-------------------------------------------------
% figure(1)
% plot(t_HH,y_HH(:,end),t_PNP,d_his,'-*')
% legend('HH model','PNP')
% title('comparison')
% grid on
% figure(2)
% plot(t_HH,y_HH(:,1:end-1),t_PNP,y(:,end-2:end),'-*')
% legend('HH model-m','HH model-n','HH model-h','PNP-m','PNP-n','PNP-h')
% title('comparison')
% grid on
% figure(3)
% plot(x,y(end,k+1:2*k),x,y(end,2*k+1:3*k),x,y(end,3*k+1:4*k))
% legend('Na','K','Cl')
% title('concentrations')
% grid on