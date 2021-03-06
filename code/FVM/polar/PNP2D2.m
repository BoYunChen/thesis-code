function dy = PNP2D2(t,y)
k = (length(y)-3)/4;
k1 = 0.5*k;
U = y(1:k);
f1 = y(k+1:2*k);
f2 = y(2*k+1:3*k);
f3 = y(3*k+1:4*k);
 m = y(end-2);
 n = y(end-1);
 h = y(end);
 
R = 1;
mid = 0.5;
dx = 0.0003;
hm = 0.005;
NA = 6.022*10^(23);
e = 1.602*10^(-19);
e0 = 8.854*10^(-12);
K = 1.38*10^(-23);
Tam = 279.45;
L = 10^(-6);
c0 = 100;
e1 = (80*e0*K*Tam)/(e*e*NA*L*L*c0);
e2 = (2*e0*K*Tam)/(e*e*NA*L*L*c0);

x1 = [mid-0.5*hm-0.5*dx:-dx:0.25*R]';
x1 = x1(end:-1:1);
x2 = [mid+0.5*hm+0.5*dx:dx:0.75*R]';
x = [x1;x2];
D = diag(1./x);


A1 = gallery('tridiag',k,1,-2,1);
A1(1,1) = -1;
A1 = e1*A1/dx/dx;
A2 = gallery('tridiag',k,-1,0,1);
A2(1,1) = -1;
A2 = e1*A2/2/dx;
A = A1 + D*A2;
dy1 = A*U + f1 + f2 - f3;

d = x(k1)*e1*(U(k1+1)-U(k1))*(1+hm)*log(x(k1+1)/x(k1))/e2/dx;
%------------------------------------------------------------------------
    g0 = 400758;
    D1 = 1.33;
    D2 = 1.96;
    D3 = 2.03;
    g11 = 1200/g0;
    g12 = 0.666/g0;
    g12 = 20*g12;
    g1 = g11*m*m*m*h + g12;
    g21 = 360/g0;
    g22 = 4.334/g0;
    g2 = g21*n^4 + g22;
    g3 = 0;
    
    J1 = D1*(0.5*(f1(2:end)+f1(1:end-1)).*(U(2:end)-U(1:end-1)) + ...
             (f1(2:end)-f1(1:end-1))).*(x(1:end-1)+x(2:end))/2/dx;
    dV1 = d - log(f1(k1)/f1(k1+1));
    J1(k1) = g1*dV1*0.5*(x1(end)+x2(1));
    J1 = [0;J1];
    df1 = (J1(2:end)-J1(1:end-1))./x(1:end-1)/dx;
    df1(k1) = df1(k1) - hm*g1*dV1/2/dx/x(k1);
    df1(k1+1) = df1(k1+1) - hm*g1*dV1/2/dx/x(k1+1);
    df1 = [df1;0];
    
    J2 = D2*(0.5*(f2(2:end)+f2(1:end-1)).*(U(2:end)-U(1:end-1)) + ...
             (f2(2:end)-f2(1:end-1))).*(x(1:end-1)+x(2:end))/2/dx;
    dV2 = d - log(f2(k1)/f2(k1+1));
    J2(k1) = g2*dV2*0.5*(x1(end)+x2(1));
    J2 = [0;J2];
    df2 = (J2(2:end)-J2(1:end-1))./x(1:end-1)/dx;
    df2(k1) = df2(k1) - hm*g2*dV2/2/dx/x(k1);
    df2(k1+1) = df2(k1+1) - hm*g2*dV2/2/dx/x(k1+1);
    df2 = [df2;0];
    
    J3 = D3*(-0.5*(f3(2:end)+f3(1:end-1)).*(U(2:end)-U(1:end-1)) + ...
             (f3(2:end)-f3(1:end-1))).*(x(1:end-1)+x(2:end))/2/dx;
    J3(k1) = 0;
    J3 = [0;J3];
    df3 = (J3(2:end)-J3(1:end-1))./x(1:end-1)/dx;
    df3 = [df3;0];
    %----------------------------------------------------------------------
V = -24.1*d+71.4;
    
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
