function uex=u_exact(x,y,k)

if (k==1)
    uex=x^2*(1-x^2)*y^2*(1-y^2);
end

if (k==2)
    uex=x*(1-x)*y*(1-y);
end
