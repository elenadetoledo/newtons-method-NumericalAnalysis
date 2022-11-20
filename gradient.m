function[grad] = gradient(x,y)

% Gradient of the himmelblaus function

dwrtx = 4*x*(x^2 + y - 11) + 2*(x + y^2 -7);
dwrty = 2*(x^2 + y - 11) + 4*y*(x + y^2 -7);

grad = [dwrtx ;dwrty];
end
