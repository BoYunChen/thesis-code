clc
clear all
close all
format long
j = 1;
while j < 18
R = 1;
x_max = 1;%0.75;
x_min = 0;%0.25;
h = 0.005;
e1 = 40;
e2 = 1;
m_m = 0.5*(R-h);
m_p = 0.5*(R+h);
N = 2*j;
[xbar,D] = legDc(N);

x1 = 0.5*(m_m-x_min)*(xbar+(m_m+x_min)/(m_m-x_min));
f1 = x1;%-exp(2*x1);%-0.01*exp(20*x1);%-cos(x1);
x2 = 0.5*(x_max-m_p)*(xbar+(x_max+m_p)/(x_max-m_p));
f2 = 1 - x2;%exp(2*(1-x2));%0.01*exp(20*(1-x2));%-cos(x2-0.01);
f = [f1(1:end-1);f2(1:end-1)];
x = [x1;x2];
%-------------u_ex--------------------------
%a1 = 0.0005;
a2 = m_m^3/3 - m_p*m_p/2 + m_p^3/3;%0.0005*(exp(20*(1-m_p))-exp(20*m_m)) + a1;
b2 = -5/36;%-0.000025 - a2;
b1 = m_p^2/4 - m_p^3/9 + a2*log(m_p) + b2 - m_m^3/9;%0.000025*(exp(20*(1-m_p))+exp(20*m_m)) + a2*m_p - a1*m_m + b2;
u_ex1 = x1.^3/9 + b1;%-0.000025*exp(20*x1) + a1*x1 + b1;
u_ex2 = x2.^2/4 - x2.^3/9 + a2*log(x2) + b2;%0.000025*exp(20*(1-x2)) + a2*x2 + b2;

jump = m_m*m_m/3;%-0.0005*exp(20*m_m) + a1;
jump = m_m*jump*e1*log(m_p/m_m)/e2;

u_ex1 = u_ex1 - jump;
u_ex = -[u_ex1;u_ex2];
%-------------------------------------------

[u,du] = SP1(40,1,f,N);
du-jump;


error(j) = max(abs(-u-u_ex));
figure(1)
plot(x,-u,'-o',x,u_ex)
set(gca,'FontSize',20)
legend('\phi_{numerical}','\phi_{exact}')
title(['N = ',num2str(N),', jump = ',num2str(-du),', error = ',num2str(error(j))])
grid on
drawnow
% figure(2);
% plot(x,u-u_ex)
% drawnow
%pause()
j = j + 1;
end
figure(2)
h=loglog(2.*[1:j-1]+1, error,'.-');
set(h,'linewidth', 2)
set(h,'MarkerSize',15)
set(gca,'FontSize',20)
set(gca,'XTick', [4 6 8 10 12 14 16])
axis([2.7 38 1e-15 0.08])
xlabel('N')
ylabel('error')
title('2D')
grid on

% u_p = -[u(1:N+1)+du;u(N+2:end)];
% figure(3)
% h = plot(x,-u,x,u_p,'-.');
% set(h, 'linewidth', 2);
% set(gca,'FontSize',30)
% set(gca, 'YTick', [],'XTick', [])
% text(0.17,-3,'d','FontSize',30)
% text(-0.01,-5.2,'0','FontSize',30)
% text(0.5,-5.2,'x_m','FontSize',30)
% text(0.99,-5.2,'1','FontSize',30)
% annotation('doublearrow',[0.25 0.25],[0.112 0.76])
% legend('\Phi','\Phi^*')
