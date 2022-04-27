function f=functionf(x,y,k)

if (k==1)
   f=-2*x^2*(1-x^2)-2*y^2*(1-y^2)+12*x^2*(1-x^2)*y^2+12*y^2*(1-y^2)*x^2;
end

if (k==2)
   f=-2*x*(x - 1) - 2*y*(y - 1);
end
