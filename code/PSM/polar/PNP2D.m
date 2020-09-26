function dy = PNP2D(t,y)
format long
k = (length(y)-1)/4;
k1 = 0.5*k;
U = y(1:k-2);
f1 = y(k-1:2*k-2);
f11 = f1(1:k1);
f12 = f1(k1+1:end);
f2 = y(2*k-1:3*k-2);
f21 = f2(1:k1);
f22 = f2(k1+1:end);
f3 = y(3*k-1:4*k-2);
f31 = f3(1:k1);
f32 = f3(k1+1:end);
 m = y(end-2);
 n = y(end-1);
 h = y(end);


NA = 6.022*10^(23);
e = 1.602*10^(-19);
e0 = 8.854*10^(-12);
K = 1.38*10^(-23);
Tam = 279.45;
L = 10^(-6);
c0 = 100;
e1 = (80*e0*K*Tam)/(e*e*NA*L*L*c0);
e2 = (2*e0*K*Tam)/(e*e*NA*L*L*c0);
rate = e1/e2;
%----------------------------------------
R = 1;
hm = 0.005;
m_m = 0.5*(R-hm);
m_p = 0.5*(R+hm);
x_min = 0;
x_max = 1;
N = k1 - 1;
[xbar,D] = legDc(N);
Dtmp = D;
Dtmp2 = 2*D/(m_m-x_min);
Dtmp3 = Dtmp2*Dtmp2;
%x1 = 0.5*m_m*(xbar+1);
x1 = 0.5*(m_m-x_min)*(xbar+(m_m+x_min)/(m_m-x_min));
%x2 = 0.5*(R-m_p)*(xbar+(R+m_p)/(R-m_p));
x2 = 0.5*(x_max-m_p)*(xbar+(x_max+m_p)/(x_max-m_p));
R1 = diag(x1);
R2 = diag(x2);
D11 = 2*R1*D/(m_m-x_min);
D12 = 2*R2*D/(x_max-m_p);
D21 = 2*diag(1./x1)*D*D11/(m_m-x_min);
D22 = 2*diag(1./x2)*D*D12/(x_max-m_p);


D = zeros(2*N,2*N);
D(1,1:N+1) = Dtmp2(1,:);
D(2:N,1:N+1) = D21(2:end-1,1:end);
D(N+1,1:N+1) = D11(end,:);
D(N+1,N+1:end) = D(N+1,N+1:end)-D12(1,1:end-1);
D(N+2:end,N+1:end) = D22(2:end-1,1:end-1);
f = f1 + f2 - f3;
f = [0;f(2:k1-1);0;f(k1+2:end-1)];
dy1 = e1*D*U + f;
% dy1tmp = Dtmp2*U(1:N+1);
% dy1tmp(1) = 0;
% dy1(1) = e1*Dtmp2(1,:)*dy1tmp + f(1);
d = Dtmp2(end,:)*U(1:N+1)*rate*x1(end)*log(x2(1)/x1(end));

U = [U(1:N);U(N+1);U(N+1:end);0];
%----------------------------------------------------------------------
    g0 = 400758;
    D1 = 1.33;
    D2 = 1.96;
    D3 = 2.03;
    g11 = 1200/g0;
    g12 = 0.65/g0;
    g1 = g11*m*m*m*h + g12;
    g21 = 360/g0;
    g22 = 4.35/g0;
    g2 = g21*n^4 + g22;
    g3 = 0;
    
    D11 = 2*R1*Dtmp/(m_m-x_min);
    D12 = 2*R2*Dtmp/(x_max-m_p);
    D21 = 2*diag(1./x1)*Dtmp/(m_m-x_min);
    D22 = 2*diag(1./x2)*Dtmp/(x_max-m_p);
    
    J11 = D1*f11.*(D11*(log(f11)+U(1:k1)));%D1*(D11*f11 + D11*U(1:k1).*f11);
    J12 = D1*f12.*(D12*(log(f12)+U(k1+1:end)));%D1*(D12*f12 + D12*U(k1+1:end).*f12);
    dV1 = d - log(f11(end)/f12(1));
    J11(1) = 0;
    J11(end) = x1(end)*g1*dV1/(1-hm);
    J12(1) = x2(1)*g1*dV1/(1+hm);
    df11 = D21*J11;
    df12 = D22*J12;
    df1 = [df11;df12];
    
%     df1tmp = D1*f11(1)*Dtmp2*(log(f11)+U(1:k1));
%     df1tmp(1) = 0;
%     df1tmp(end) = g1*dV1;
%     df1(1) = 2*Dtmp2(1,:)*df1tmp;

    df1(1) = Dtmp2(1,:)*f11;
    df1(end) = 0;
    %df1(1) = 0;
    
    J21 = D2*f21.*(D11*(log(f21)+U(1:k1)));%D2*(D11*f21 + D11*U(1:k1).*f21);
    J22 = D2*f22.*(D12*(log(f22)+U(k1+1:end)));%D2*(D12*f22 + D12*U(k1+1:end).*f22);
    dV2 = d - log(f21(end)/f22(1));
    J21(1) = 0;
    J21(end) = x1(end)*g2*dV2/(1-hm);
    J22(1) = x2(1)*g2*dV2/(1+hm);
    df21 = D21*J21;
    df22 = D22*J22;
    df2 = [df21;df22];
    
%     df2tmp = D2*f21(1)*Dtmp2*(log(f21)+U(1:k1));
%     df2tmp(1) = 0;
%     df2tmp(end) = g2*dV2;
%     df2(1) = 2*Dtmp2(1,:)*df2tmp;

    df2(1) = Dtmp2(1,:)*f21;
    df2(end) = 0;
    %df2(1) = 0;
    
    J31 = D3*f31.*(D11*(log(f31)-U(1:k1)));%D3*(D11*f31 - D11*U(1:k1).*f31);
    J32 = D3*f32.*(D12*(log(f32)-U(k1+1:end)));%D3*(D12*f32 - D12*U(k1+1:end).*f32);
    dV3 = d - log(f31(end)/f32(1));
    J31(1) = 0;
    J31(end) = x1(end)*g3*dV3;
    J32(1) = x2(1)*g3*dV3;
    df31 = D21*J31;
    df32 = D22*J32;
    df3 = [df31;df32];
    
%     df3tmp = D3*f31(1)*Dtmp2*(log(f31)-U(1:k1));
%     df3tmp(1) = 0;
%     df3tmp(end) = g3*dV3;
%     df3(1) = 2*Dtmp2(1,:)*df3tmp;

    df3(1) = Dtmp2(1,:)*f31;
    df3(end) = 0;
    %df3(1) = 0;
 
    %----------------------------------------------------------------------
V = -24.07*d+71.4;
    
dim = 10^(0);
    
a_m = 0.1*(V-25)/(1-exp((25-V)/10))*dim;
b_m = 4*exp((-V)/18)*dim;
a_h = 0.07*exp((-V)/20)*dim;
b_h = 1/(1+exp((30-V)/10))*dim;
a_n = 0.01*(V-10)/(1-exp(-(V-10)/10))*dim;
b_n = 0.125*exp((-V)/80)*dim;
    
    
m_t = a_m*(1-m) - b_m*m;
n_t = a_n*(1-n) - b_n*n;
h_t = a_h*(1-h) - b_h*h;
    
dy = [dy1;df1;df2;df3;m_t;n_t;h_t];
end