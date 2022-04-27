function uex=u_exact(x,y,k)

if (k==1)
    uex=cos(2*pi*x)*sin(2*pi*y^2);
end

if (k==2)
    uex=y^2*(1-y);
end
