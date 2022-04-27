function f=functionf(x,y,k)

if (k==1)
    f=4*pi*cos(2*pi*x)*(pi*sin(2*pi*y^2)*(1+4*y^2)-cos(2*pi*y^2));
end

if (k==2)
    f=6*y-2;
end
