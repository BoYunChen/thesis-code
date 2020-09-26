clear all
close all
clc
format long

N = 70;
hm = 0.005;
R = 1;
mid = 0.5;
T1 = 100;
T2 = 25;

[x,f1,f2,f3] = createf(N,hm,R,mid);
U = zeros(length(f1)-2,1);

n = 0.3177;
m = 0.05293;
h = 0.5961;

y = [U;f1;f2;f3;m;n;h];
k = (length(y)-1)/4;
A = ones(1,length(y));
%A(1:k-2) = 0*A(1:k-2);

A(1:k-1) = 0*A(1:k-1);
A(2*k-1) = 0;
A(3*k-1) = 0;

A = sparse(diag(A));
options = odeset('Mass', A);
[t, y]= ode15s(@PNP1D, [0 T1], y,options);
[t1, y1] = ode15s(@PNP1D2,[0 0.1],y(end,:),options);
[t2, y2] = ode15s(@PNP1D,[0.1 T2],y1(end,:),options);

t = [t1(1:end-1);t2];
y = [y1(1:end-1,:);y2];

U = y(:,1:k-2);
d = fVm(y(:,1:0.5*k)',N);
U = [U(:,1:0.5*k) U(:,0.5*k:end) 0*U(:,end)];


c1 = y(:,k-1:2*k-2);
c2 = y(:,2*k-1:3*k-2);
c3 = y(:,3*k-1:4*k-2);
m = y(:,4*k-1);
n = y(:,4*k);
h = y(:,4*k+1);

t_HH = csvread('t_HH.csv');
y_HH = csvread('y_HH.csv');

j = 1;
while j < length(t)+1
    figure(1)
    h1=plot(x,c1(j,:),x,c2(j,:),x,c3(j,:));
    set(h1, 'linewidth', 2);
    title(['t=',num2str(t(j)),'/',num2str(T2)])
    set(gca,'FontSize',25)
    grid on
    drawnow
    
    figure(2)
    h2=plot(x,c1(j,:)+c2(j,:)-c3(j,:));
    set(h2, 'linewidth', 2);
    set(gca,'FontSize',25)
    grid on
    drawnow
    
    U(j,1:0.5*k) = U(j,1:0.5*k) - d(j);
    figure(3)
    h3 = plot(x,U(j,:));
    set(h3, 'linewidth', 2);
    set(gca,'FontSize',25)
    grid on
    pause(0.1)
    j = j + 1;   
end
y = [U y(:,k-1:end)];
dlmwrite('t_1dsp70.csv',t,'precision',10)
dlmwrite('y_1dsp70.csv',y,'precision',10)
dlmwrite('x_1dsp70.csv',x,'precision',10)

figure(4)
h4 = plot(t,-24.1*d,'-o',t_HH,y_HH(:,end));
set(h4, 'linewidth', 2);
xlabel('t(ms)')
ylabel('Membrane potential(mV)')
set(gca,'FontSize',25)
legend('PNP','HH')
grid on

% figure(4)
% plot(t,m,t,n,t,h)
% legend('m','n','h')
% grid on

% y0 = y(end,:);
% 
% [t1, y1] = ode15s(@PNP1D2,[0 0.1],y0,options);
% d_his1 = 40*hm*(y1(:,0.5*k+1)-y1(:,0.5*k))/dx;
% 
% [t2, y2] = ode15s(@PNP1D,[0.1:20],y1(end,:),options);
% d_his2 = 40*hm*(y2(:,0.5*k+1)-y2(:,0.5*k))/dx;
% j = 1;
% while j < length(t2)+1
%     y2(j,1:0.5*k) = y2(j,1:0.5*k) - d_his2(j);
%     figure(1)
%     subplot(1,2,1),plot(x,24.1*y2(j,1:k),'-.')
%     grid on
%     title(['t=',num2str(t(j)),'/',num2str(T)])
%     %-------------------------------------------
%     subplot(1,2,2),plot(x,y2(j,k+1:2*k),x,y2(j,2*k+1:3*k),x,y2(j,3*k+1:4*k))
%     legend('Na^+','K^+','Cl^-')
%     grid on
%     title(['t=',num2str(t(j)),'/',num2str(T)])
%     %------------------------------------------------
%     pause()
%     j = j+1;
% end

%---------------------------------------------------------------

%net = y(:,k+1:2*k) + y(:,2*k+1:3*k) - y(:,3*k+1:4*k);
%d_his = 40*hm*(y(:,0.5*k+1)-y(:,0.5*k))/dx;
%y(:,1:0.5*k) = y(:,1:0.5*k) - d_his;
%plot(x,24.1*y(end,1:k))

% j = 1;
% while j < length(t)+1
%     figure(1)
%     subplot(1,2,1),plot(x,y(j,1:k) - d_his(j,:))
%     grid on
%     title(['t=',num2str(t(j)),'/',num2str(T)])
%     %-------------------------------------------
%     subplot(1,2,2),plot(x,y(j,k+1:2*k),x,y(j,2*k+1:3*k),x,y(j,3*k+1:4*k))
%     legend('Na^+','K^+','Cl^-')
%     grid on
%     title(['t=',num2str(t(j)),'/',num2str(T)])
%     %------------------------------------------------
% %     subplot(2,2,3),plot(x,net(j,:),'-o')
% %     grid on
% %     title(['t=',num2str(t(j)),'/',num2str(T)])
%     %--------------------------------------------------
%     pause()
%     j = j+1;
% end
% c_his = [y(:,1.5*k+1)./y(:,1.5*k) y(:,2.5*k+1)./y(:,2.5*k)]';
% c_his = 24.1*log(c_his);
% figure(2)
% subplot(2,2,1),plot(t,-d_his,'-o')
% title('membrane potential')
% grid on
% subplot(2,2,2),plot(t,c_his(1,:),'-o',t,c_his(2,:),'-*')
% legend('Na^+','K^+')
% title('Nernst potential')
% grid on
% subplot(2,2,3),plot(t,y(:,end-2),t,y(:,end-1),t,y(:,end))
% legend('m','n','h')
% title('gatting variable')
% grid on
%saveas(4,'Vm_10.m')