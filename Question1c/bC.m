function g=bC(x,y,k)

if k==1
    g=cos(2*pi*x)*sin(2*pi*y^2);
end

if k==2
    g=y^2*(1-y);
end