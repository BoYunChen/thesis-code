function [x,f1,f2,f3] = createf(dx,hm,R,m)
x1 = [m-0.5*hm-0.5*dx:-dx:0.25*R]';
x1 = x1(end:-1:1);
k1 = length(x1);
x2 = [m+0.5*hm+0.5*dx:dx:0.75*R]';
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