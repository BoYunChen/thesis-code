function V = fVm(U,t,N)
R = 1;
hm = 0.005;
m_m = 0.5*(R-hm);
m_p = 0.5*(R+hm);
[xbar,D] = legDc(N);
D11 = 2*D/m_m;
D12 = 2*D/(R-m_p);
D21 = D11*D11;
D22 = D12*D12;
rate = 40;
D = zeros(2*N,2*N);
D(1,1:N+1) = D11(1,:);
D(2:N,1:N+1) = D21(2:end-1,1:end);
D(N+1,1:N+1) = D11(end,:);
D(N+1,N+1:end) = D(N+1,N+1:end)-D12(1,1:end-1);
D(N+2:end,N+1:end) = D22(2:end-1,1:end-1);
hm = lsf(0,25,hm,0.5,0.05,t);
%hm = usf(0,25,hm,2,0.1,t);
% size(D11(end,:))
% size(U)
% size(hm)
V = D11(end,:)*U*rate.*hm';
end