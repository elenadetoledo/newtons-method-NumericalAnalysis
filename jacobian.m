function [jaco] = jacobian( x, y)

% Jacobian matrix from the himmelblaus function

Comp1 = 12*(x^2)+4*y-42;
Comp2 = 4*x+4*y;
Comp3 = 12*(y^2)+4*x-26;

jaco = [ Comp1, Comp2 ; Comp2, Comp3];

end

