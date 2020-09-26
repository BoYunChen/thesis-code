function y = lsf(tmin,tmax,hm,scale, epsilon,t)
y = (1-scale)*hm./(1+exp((t-0.5*(tmax-tmin))./epsilon));
y = y + scale*hm;
end