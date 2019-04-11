% Multiply the 4-by-4 matrix M by the vector [x1, x2, x3, 1].
function [y1,y2,y3,y4] = coords(M, x1, x2, x3)

y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
y4 = M(4,1)*x1 + M(4,2)*x2 + M(4,3)*x3 + M(4,4);
