function v1 = vanishing_point(xo1, xf1, xo2, xf2)
    x = cross(xo1,xf1); y = cross(xo2,xf2);
    v1 = cross(x/x(end),y/y(end)); 
end