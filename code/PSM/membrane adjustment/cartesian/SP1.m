function [U,du] = SP1(N,f)
R = 1;
h = 0.005;
alpha = 0; %u'(0)
beta = 0;  %u(R)
rate = 40;
m_m = 0.5*(R-h);
m_p = 0.5*(R+h);
[xbar,D] = legDc(N);
f1 = f(1:N);
f2 = f(N+1:2*N);
f = [f1;f2];
D11 = 2*D/m_m;
D12 = 2*D/(R-m_p);
D21 = D11*D11;
D22 = D12*D12;



D = zeros(2*N,2*N);
D(1,1:N+1) = D11(1,:);
D(2:N,1:N+1) = D21(2:end-1,1:end);
D(N+1,1:N+1) = D11(end,:);
D(N+1,N+1:end) = D(N+1,N+1:end)-D12(1,1:end-1);
D(N+2:end,N+1:end) = D22(2:end-1,1:end-1);

F=f;
F(1) = alpha;
F(N+1:end) = F(N+1:end) - [-D12(1,end);D22(2:end-1,end)]*beta;
F(N+1) = 0;
u = D\F;
u = [u;beta];
du = D11(end,:)*u(1:N+1)*rate*h;
U = [u(1:N+1)-du;u(N+1);u(N+2:end)];
%U = -[u(1:N+1)+du;u(N+2:end)];
end