clear all
close all
clc
format long

N = 70;
hm = 0.005;
R = 1;
x_min = 0;
x_max = 1;
mid = 0.5;
T1 = 30;
T2 = 25;

[x,f1,f2,f3] = createf(N,hm,x_min,x_max,R);
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
options2 = odeset('Mass', A,'AbsTol',5e-7);
[t, y]= ode15s(@PNP2D, [0 T1], y,options);
[t1, y1] = ode15s(@PNP2D2,[0 0.1],y(end,:),options);
[t2, y2] = ode15s(@PNP2D,[0.1 T2],y1(end,:),options);

t = [t1(1:end-1);t2];
y = [y1(1:end-1,:);y2];

U = y(:,1:k-2);
d = fVm(y(:,1:0.5*k)',N,x_min,x_max);
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
%aviobj = VideoWriter('COM60.avi','Uncompressed AVI');
%aviobj.FrameRate = 5;
%open(aviobj)
 while j < length(t)+1
    figure(j)
    subplot(2,2,1), h11=plot(t,-24.07*d,t(j),-24.07*d(j),'-or');
    axis([0 T2 -82 41])
    xlabel('t(ms)')
    ylabel('membrane potential(mV)')
    set(h11, 'linewidth', 2);
    set(gca,'FontSize',20)
    
    subplot(2,2,2), h12 = plot(x,100*c1(j,:),x,100*c2(j,:),x,100*c3(j,:));
    xlabel('x(\mu m)')
    ylabel('concentration(mM)')
    legend('Na','K','Cl','Location','Best');
    set(h12, 'linewidth', 1.5);
    set(gca,'FontSize',20)
    
    U(j,1:0.5*k) = U(j,1:0.5*k) - d(j);
    subplot(2,2,3), h21 = plot(x,24.07*U(j,:));
    xlabel('x(\mu m)')
    ylabel('electric potential(mV)')
    axis([0 1 -87 45])
    set(h21, 'linewidth', 2);
    set(gca,'FontSize',20)
    
    subplot(2,2,4), h22=plot(x,100*(c1(j,:)+c2(j,:)-c3(j,:)));
    xlabel('x(\mu m)')
    ylabel('net concentration(mM)')
    axis([0 1 -3.5 3.5])
    set(h22, 'linewidth', 2);
    set(gca,'FontSize',20)
    
    set(gcf,'outerposition',get(0,'screensize'));
    
    
    frames1(j)=getframe(gcf);
    drawnow
    close all
%     figure(j)
%     subplot(2,1,1),h11=plot(t,-24.1*d,t(j),-24.1*d(j),'-or');
%     set(h11, 'linewidth', 2);
%     subplot(2,1,2),h12=plot(x,100*c1(j,:),x,100*c2(j,:),x,100*c3(j,:));
%     axis([0.45 0.55 0 150])
%     set(h12, 'linewidth', 2);
%     %title(['t=',num2str(t(j)),'/',num2str(T2)])
%     set(gca,'FontSize',25)
%     grid on
%     frames1(j)=getframe(gcf);
%     close all
    
    
%     drawnow
%     
%     figure(j)
%     subplot(2,1,1),h21=plot(t,-24.1*d,t(j),-24.1*d(j),'-or');
%     set(h21, 'linewidth', 2);
%     subplot(2,1,2),h22=plot(x,100*(c1(j,:)+c2(j,:)-c3(j,:)));
%     set(h22, 'linewidth', 2);
%     axis([0.45 0.55 -3.3 3.3])
%     set(gca,'FontSize',25)
%     xlabel('x(\mu m)')
%     ylabel('concentrations(mM)')
%     dim = [0.125 0.49 0.05 0.2];
%     annotation('rectangle',dim,'Color','red')
%     title(['t=',num2str(t(j)),'/',num2str(T2)])
%     grid on
%     frames2(j)=getframe(gcf);
%     close all
%     
%     drawnow
    
%      U(j,1:0.5*k) = U(j,1:0.5*k) - d(j);
%     figure(3)
%     subplot(2,1,1),h31=plot(t,-24.1*d,t(j),-24.1*d(j),'-or');
%     set(h31, 'linewidth', 2);
%     subplot(2,1,2),h32 = plot(x,24.1*U(j,:));
%     set(h32, 'linewidth', 2);
%     axis([0.45 0.55 -87 45])
%     %set(h3, 'linewidth', 2);
%     set(gca,'FontSize',25)
%     %title(['t=',num2str(t(j)),'/',num2str(T2)])
%     grid on
%     frames3(j)=getframe(gcf);
%     close all
     pause(0.1)
     j = j + 1;
     
 end
dt = 0.1;
genGIF(length(t),frames1,dt,'concnetration70f.gif')
% writeVideo(aviobj,frames1)
% close(aviobj)
% genGIF(length(t),frames2,dt,'net_charge70_1.gif')
% genGIF(length(t),frames3,dt,'potential70_1.gif')

 y = [U y(:,k-1:end)];
 dlmwrite('t_2dsp70.csv',t,'precision',10)
 dlmwrite('y_2dsp70.csv',y,'precision',10)
 dlmwrite('x_2dsp70.csv',x,'precision',10)

%-----------------------------------------------------------------

% figure(1)
% h1 = plot(x,24.07*[U(51,1:0.5*k)-d(51) U(51,0.5*k+1:end)]);
% axis([0 1 0 42])
% set(h1, 'linewidth', 2);
% xlabel('x(\mu m)')
% ylabel('electric potential(mV)')
% set(gca,'FontSize',20)

% figure(1)
% h1 = plot(x,24.07*[U(51,1:0.5*k)-d(51) U(51,0.5*k+1:end)]);
% axis([0.48 0.4978 40.4 40.5])
% set(h1, 'linewidth', 2);
% set(gca,'FontSize',25)

% figure(1)
% h1 = plot(x,24.07*[U(51,1:0.5*k)-d(51) U(51,0.5*k+1:end)]);
% axis([0.5022 0.52 -0.01 0.3])
% set(h1, 'linewidth', 2);
% set(gca,'FontSize',25)


% figure(2)
% h2 = plot(x,100*c1(51,:),x,100*c2(51,:),x,100*c3(51,:));
% set(h2, 'linewidth', 2);
% axis([0.4 0.6 0 140])
% legend('Na','K','Cl','Location','Best')
% xlabel('x(\mu m)')
% ylabel('concentration(mM)')
% set(gca,'FontSize',20)
% 
% figure(3)
% h3 = plot(x,100*(c1(51,:)+c2(51,:)-c3(51,:)));
% set(h3, 'linewidth', 2);
% axis([0.4 0.6 -2 2])
% xlabel('x(\mu m)')
% ylabel('net concentration(mM)')
% set(gca,'FontSize',20)
