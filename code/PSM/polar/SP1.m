function [U,du] = SP1(e1,e2,f,N)
R = 1;
h = 0.005;
alpha = 0; %u'(0)
beta = 0;  %u(R)
rate = e1/e2;
x_max = 1;%0.75;
x_min = 0;%0.25;
m_m = 0.5*(R-h);
m_p = 0.5*(R+h);
[xbar,D] = legDc(N);
Dtmp = 2*D/(m_m-x_min);
Dtmp2 = Dtmp*Dtmp;
f1 = f(1:N);
f2 = f(N+1:2*N);
f = [f1;f2];
x1 = 0.5*(m_m-x_min)*(xbar+(m_m+x_min)/(m_m-x_min));
x2 = 0.5*(x_max-m_p)*(xbar+(x_max+m_p)/(x_max-m_p));
R1 = diag(x1);
R2 = diag(x2);
D11 = 2*R1*D/(m_m-x_min);
D12 = 2*R2*D/(x_max-m_p);
D21 = 2*diag(1./x1)*D*D11/(m_m-x_min);
D22 = 2*diag(1./x2)*D*D12/(x_max-m_p);


D = zeros(2*N,2*N);
D(1,1:N+1) = Dtmp(1,:);
D(2:N,1:N+1) = D21(2:end-1,1:end);
D(N+1,1:N+1) = D11(end,:);
D(N+1,N+1:end) = D(N+1,N+1:end)-D12(1,1:end-1);
D(N+2:end,N+1:end) = D22(2:end-1,1:end-1);

F=f;
F(1) = 0;%0.5*F(1);
F(N+1) = 0;
u = D\F;
u = [u;beta];
du = Dtmp(end,:)*u(1:N+1)*rate*x1(end)*log(x2(1)/x1(end));%du = D11(end,:)*u(1:N+1)*rate*h;
U = [u(1:N+1)-du;u(N+1:end)];
%U = -[u(1:N+1)+du;u(N+2:end)];
end