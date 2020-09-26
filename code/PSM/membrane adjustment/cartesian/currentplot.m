function currentplot(t,d,c1,c2,m,n,h)
    k = length(c1);
    g0 = 400758;
    g11 = 1200/g0;
    g12 = 0.65/g0;
    g21 = 360/g0;
    g22 = 4.35/g0;
    Na_main = -g11*(d'-log(c1(:,0.5*k)./c1(:,0.5*k+1)));
    Na_main = m.^3.*h.*Na_main;
    Na_leak = -g12*(d'-log(c1(:,0.5*k)./c1(:,0.5*k+1)));
    K_main = -g21*(d'-log(c2(:,0.5*k)./c2(:,0.5*k+1)));
    K_main = n.^4.*K_main;
    K_leak = -g22*(d'-log(c2(:,0.5*k)./c2(:,0.5*k+1)));
    
    figure(1)
    h1 = plot(t,Na_main,t,K_main);
    set(h1, 'linewidth', 2);
    title('main current')
    legend('Na','K')
    set(gca,'FontSize',25)
    grid on
    figure(2)
    h2 = plot(t,Na_leak,t,K_leak);
    set(h2, 'linewidth', 2);
    title('leak current')
    legend('Na','K')
    set(gca,'FontSize',25)
    grid on
    drawnow
end