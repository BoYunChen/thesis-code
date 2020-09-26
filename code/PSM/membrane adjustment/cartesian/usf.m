function y = usf(tmin,tmax,hm,scale, epsilon,t)
y = (scale-1)*hm./(1+exp(-(t-0.5*(tmax-tmin))./epsilon));
y = y + hm;
end