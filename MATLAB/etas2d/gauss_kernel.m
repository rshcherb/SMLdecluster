function fVal = gauss_kernel(r,bw)
%
%   Computes a spherically symmetric Gaussian kernel value 
%
%   r  - distance from the centre
%   bw - bandwidth
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 26 May 2021
%   ...
%   version: 1.0.0, 26 May 2021
%
    bw2 = 1.0/(2.0*bw^2);
    fVal = bw2/pi*exp(-r^2.*bw2);             % Gaussian kernel function normlaized to 1 
    %fVal = 1.0/(2.0*pi*bw)*exp(-r^2.*bw2);   % Gaussian kernel function from Zhuang, EPS 2011, Eq. (10) 
end

