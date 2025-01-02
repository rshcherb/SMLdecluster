function [x, y] = dr_rad_inverse(Px,Py,d,q)
%
%   Computes the random distance from the parent event acording to Zhuang et al. (2005), Model 3 and 5
%   normalization: (q-1)/(\pi*d^2)*(1 + r^2/d^2)^(-q)
%   Zhuang and Tuati (CORSSA V, 2012, p.8)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 April 2021
%   ...
%   version: 1.0.0, 10 May 2021
%
    U = 0.0;
    while U == 0.0
        U = rand();
    end
    
    r = d*sqrt(U^(1/(1-q)) - 1);
    
    phi = 2*pi*rand();
    x = Px + r*cos(phi);
    y = Py + r*sin(phi);
end
