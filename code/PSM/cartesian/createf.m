function [x,f1,f2,f3] = createf(N,h,R,m)

m_m = 0.5*(R-h);
m_p = 0.5*(R+h);
[xbar,D] = legDc(N);
x1 = 0.5*m_m*(xbar+1);
x2 = 0.5*(R-m_p)*(xbar+(R+m_p)/(R-m_p));
x = [x1;x2];

f11 = 0*x1 + 0.12;
f12 = 0*x2 + 1;%1
f21 = 0*x1 + 1.25;
f22 = 0*x2 + 0.04;%0.04
f31 = 0*x1 + 1.37;
f32 = 0*x2 + 1.04;%1.04

f1 = [f11;f12];
f2 = [f21;f22];
f3 = [f31;f32];
end