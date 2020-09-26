clear all
close all

hm = 0.0005;
scale = 0.5;
epsilon = 0.05;
dx = 0.001;
x = [0:dx:25];
sf = (1-scale)*hm./(1+exp((x-12.5)./epsilon));
sf = sf + scale*hm;
plot(x,sf,'-o')