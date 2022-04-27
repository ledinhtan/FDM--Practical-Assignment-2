function g=bC(x,y,k,h)

if k==1
    if h==1
        g=sin(2*pi*(x-1));
    elseif h==2
        g=-10*y^2;
    elseif h==3
        g=sin(2*pi*(x-1))-10;
    else
        g=2*pi;
    end
end

if k==2
    if h==1
        g=x^2;
    elseif h==2
        g=y^2;
    elseif h==3
        g=x^2+1;
    else
        g=2;
    end
end

if k==3
    if h==1
        g=4/9*(x-1/3)^2;
    elseif h==2
        g=1/9*(y-2/3)^2;
    elseif h==3
        g=1/9*(x-1/3)^2;
    else
        g=4/3*(y-2/3)^2;
    end
end




