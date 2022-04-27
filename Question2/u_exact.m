function uex=u_exact(x,y,k)

if (k==1)
    uex=sin(2*pi*(x - 1)) - 10*y^2;
end

if (k==2)
    uex=x^2+y^2;
end

if (k==3)
    uex=(x-1/3)^2*(y-2/3)^2;
end
