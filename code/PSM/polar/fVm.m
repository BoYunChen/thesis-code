function V = fVm(U,N,x_min,x_max)
R = 1;
hm = 0.005;
rate = 40;
m_m = 0.5*(R-hm);
m_p = 0.5*(R+hm);
[xbar,D] = legDc(N);
Dtmp = D;
Dtmp2 = 2*Dtmp/(m_m-x_min);
x1 = 0.5*(m_m-x_min)*(xbar+(m_m+x_min)/(m_m-x_min));
x2 = 0.5*(x_max-m_p)*(xbar+(x_max+m_p)/(x_max-m_p));

V = Dtmp2(end,:)*U*rate*x1(end)*log(x2(1)/x1(end));
end