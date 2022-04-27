function uex=u_exact(x,y,k)

if (k==1)
    uex=cos(2*pi*x)*cos(5*pi*y);
end

if (k==2)
    uex=cos(3*pi*x)*cos(5*pi*y);
end

if (k==3)
    uex = (x^2/2 - x^3/3)*(y^2/2 - y^3/3)-1/144;
end

