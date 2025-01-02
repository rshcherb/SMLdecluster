function [vCx, vCy] = triangle_point(Ax,Ay,Bx,By,a,b,c)
%
%   Computes the coordinates of the third point of a tringle if two points
%   (Ax,Ay) (Bx,By) and the sides a, b, c are known
%   http://paulbourke.net/geometry/circlesphere/
%   https://math.stackexchange.com/questions/187107/calculate-coordinates-of-3rd-point-vertex-of-a-scalene-triangle-if-angles-and?rq=1
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 29 May 2021
%   ...
%   version: 1.0.0, 29 May 2021
%
%           (Cx,Cy)
%          /|
%        b/ |
%        /  |
%(Ax,Ay) \  |a
%        c\ |
%          \|
%           (Bx,By)
%            

    ha = (b^2 - a^2 + c^2)/(2*c);

    hp = sqrt(b^2 - ha^2);
    hm = -sqrt(b^2 - ha^2);
    
    x2 = Ax + ha*(Bx - Ax)/c;
    y2 = Ay + ha*(By - Ay)/c;
    
    vCx(1) = x2 + hm*(By - Ay)/c;
    vCy(1) = y2 + hp*(Bx - Ax)/c;

    vCx(2) = x2 + hp*(By - Ay)/c;
    vCy(2) = y2 + hm*(Bx - Ax)/c;
end
