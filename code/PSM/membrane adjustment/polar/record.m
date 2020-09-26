clear all
close all

format long 

x1 = csvread('x_1dsp70.csv');
t1 = csvread('t_1dsp70.csv');
y1 = csvread('y_1dsp70.csv');

x2 = csvread('x_2dsp70.csv');
t2 = csvread('t_2dsp70.csv');
y2 = csvread('y_2dsp70.csv');


t_HH = csvread('t_HH.csv');
y_HH = csvread('y_HH.csv');

T1 = t1(end);
T2 = t2(end);
k1 = length(x1);
k2 = length(x2);

d_his1 = y1(:,0.5*k1+1) - y1(:,0.5*k1);
d_his2 = y2(:,0.5*k2+1) - y2(:,0.5*k2);

re_Na1 = log(y1(:,1.5*k1+1)./y1(:,1.5*k1));
re_K1 =  log(y1(:,2.5*k1+1)./y1(:,2.5*k1));

re_Na2 = log(y2(:,1.5*k2+1)./y2(:,1.5*k2));
re_K2 =  log(y2(:,2.5*k2+1)./y2(:,2.5*k2));

j = 1;

figure(1)
h1=plot(t1,-24.1*d_his1,'-ob',t2,-24.1*d_his2,'-or',t_HH,y_HH(:,end));
set(h1, 'linewidth', 2);
xlabel('t(ms)')
ylabel('Membrane potential(mV)')
set(gca,'FontSize',25)
legend('PNP-Cartesian','PNP-polar','HH','Location','Best')
grid on

figure(4)
h4=plot(t1,24.1*re_Na1,t2,24.1*re_Na2);
title('Reversal potential for Na^+')
set(h4, 'linewidth', 2);
xlabel('t(ms)')
ylabel('Reversal potential(mV)')
%axis([0 20 -90 60])
set(gca,'FontSize',25)
legend('Cartesian','polar','Location','Best')
grid on

figure(5)
h5=plot(t1,24.1*re_K1,t2,24.1*re_K2);
title('Reversal potential for K^+')
set(h5, 'linewidth', 2);
xlabel('t(ms)')
ylabel('Reversal potential(mV)')
%axis([0 20 -90 60])
set(gca,'FontSize',25)
legend('Cartesian','polar','Location','Best')
grid on



% figure(1)
% h11=plot(t2,-24.1*d_his2);
% hold on
% h12=plot(t2(53),-24.1*d_his2(53),'o');
% set(h11, 'linewidth', 2);
% set(h12, 'linewidth', 4);
% xlabel('t(ms)')
% ylabel('Membrane potential(mV)')
% set(gca,'FontSize',25)
% 
% figure(2)
% h2 = plot(x2,24.1*y2(53,1:k2));
% set(h2,'linewidth',2)
% xlabel('x(\mu m)')
% ylabel('electric potential(mV)')
% set(gca,'FontSize',25)
% 
% figure(3)
% h3 = plot(x2,100*y2(53,k2+1:2*k2)...
%          ,x2,100*y2(53,2*k1+1:3*k2)...
%          ,x2,100*y2(53,3*k2+1:4*k2));
% legend('Na','K','Cl')
% set(h3,'linewidth',2)
% xlabel('x(\mu m)')
% ylabel('concentrations(mM)')
% set(gca,'FontSize',25)
% 
% figure(4)
% h4 = plot(x2,100*y2(53,k2+1:2*k2)...
%          +100*y2(53,2*k1+1:3*k2)...
%          -100*y2(53,3*k2+1:4*k2));
% set(h4,'linewidth',2)
% xlabel('x(\mu m)')
% ylabel('net concentrations(mM)')
% set(gca,'FontSize',25)

% j = 1;
% while j < length(t)+1
% %     figure(1)
% %     subplot(1,2,1),plot(x,24.1*y(j,1:k),'-o')
% %     grid on
% %     title(['t=',num2str(t(j)),'/',num2str(T)])
% %     %-------------------------------------------
% %     subplot(1,2,2),plot(x,y(j,k+1:2*k),x,y(j,2*k+1:3*k),x,y(j,3*k+1:4*k))
% %     legend('Na^+','K^+','Cl^-')
% %     grid on
% %     title(['t=',num2str(t(j)),'/',num2str(T)])
% %     drawnow
%     j = j + 1;
% end
% figure(1)
% h1=plot(x1,24.1*y1(1,1:k1),x2,24.1*y2(1,1:k2),'o');
% legend('1D','2D')
% set(h1, 'linewidth', 2);
% xlabel('x(\mu m)')
% ylabel('Electric potential(mV)')
% set(gca,'FontSize',25)
% %grid on
% 
% figure(2)
% h2=plot(t1,-24.1*d_his1,'-ob',t2,-24.1*d_his2,'-or',t_HH,y_HH(:,end));
% set(h2, 'linewidth', 2);
% xlabel('t(ms)')
% ylabel('Membrane potential(mV)')
% set(gca,'FontSize',25)
% legend('PNP-1D','PNP-2D','HH')
% grid on
% figure(3)
% h3=plot(x1,100*(y1(end,k1+1:2*k1)+y1(end,2*k1+1:3*k1)-y1(end,3*k1+1:4*k1)),...
%         x2,100*(y2(end,k2+1:2*k2)+y2(end,2*k2+1:3*k2)-y2(end,3*k2+1:4*k2)),'o');
% legend('1D','2D')
% axis([0.4 0.6 -4 3])
% set(h3, 'linewidth', 2);
% set(gca,'FontSize',20)
% xlabel('x(\mu m)')
% ylabel('concentrations(mM)')
% grid on
% figure(4)
% h4=plot(t1,24.1*re_Na1,t2,24.1*re_Na2);
% title('Reversal potential for Na^+')
% set(h4, 'linewidth', 2);
% xlabel('t(ms)')
% ylabel('Reversal potential(mV)')
% %axis([0 20 -90 60])
% set(gca,'FontSize',25)
% legend('1D','2D')
% grid on
% 
% figure(5)
% h5=plot(t1,24.1*re_K1,t2,24.1*re_K2);
% title('Reversal potential for K^+')
% set(h5, 'linewidth', 2);
% xlabel('t(ms)')
% ylabel('Reversal potential(mV)')
% %axis([0 20 -90 60])
% set(gca,'FontSize',25)
% legend('1D','2D')
% grid on