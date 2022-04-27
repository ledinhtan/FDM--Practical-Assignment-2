function uex=u_exact(x,y,k)

if (k==1)
    uex=x*(1-x);
end

if (k==2)
    uex=x^2*(1-x);
end
