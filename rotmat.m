% Compose the 4x4 matrix corresponding to a rotation about the x-, y- and
% z-axes (specified in radians).
function mat = rotmat(ang)

if numel(ang) < 3
    ang(end+1:3) = 0;
end

xrot = [ 1         0            0            0;
         0         cos(ang(1)) -sin(ang(1))  0;
         0         sin(ang(1))  cos(ang(1))  0;
         0         0            0            1];
  
yrot = [ cos(ang(2))  0         sin(ang(2))  0;
         0            1         0            0;
        -sin(ang(2))  0         cos(ang(2))  0;
         0            0         0            1];

zrot = [ cos(ang(3)) -sin(ang(3))  0         0;
         sin(ang(3))  cos(ang(3))  0         0;
         0            0            1         0;
         0            0            0         1];

mat = xrot * yrot * zrot;
