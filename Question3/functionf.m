function f=functionf(x,y,k)

if (k==1)
    f=29*pi^2*cos(2*pi*x)*cos(5*pi*y);
end

if (k==2)
    f=34*pi^2*cos(3*pi*x)*cos(5*pi*y);
end

if (k==3)
    f = (2*x -1)*(y^2/2-y^3/3) + (2*y -1)*(x^2/2-x^3/3);
end



