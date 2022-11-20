function[Next_x] = newtonraphson( x, y);

% We need to compute the arguments of our function Xn+1 = Xn - jacobian^(-1)*gradient

product = inv(jacobian(x,y))*gradient(x,y);
PointX = x - product(1);
PointY = y - product(2);
Next_x = [ PointX , PointY ];
end