function [mU, vFlist] = faults_intensity(vX,vY,nF,vFpar,fH,varargin)
%
%   Generates a background earthquake rate u(x,y) on a rectangular grid with presence of faults
%
%   vX, vY - the size of the grid
%   nF     - number of faults to generate
%   fH     - the Hurst exponent: (0, 1)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 27 May 2021
%   ...
%   version: 1.0.0, 17 June 2021
%
    bNormFactor = false;
    for k = 1:length(varargin)
        if strcmp('NormFactor',varargin{k})
            bNormFactor = true;
            fNormFactor = varargin{k+1};
        end
    end
    q   = vFpar(1);    % parameter that controls the distrobution of fault lengths
    d   = vFpar(2);    % second parameter for the distribution of fault lengths
    fP  = vFpar(3);    % mixing parameter for combining fractal and fault intensities
    fI0 = vFpar(4);    % lower level below which the intensity is set to zero.
    nX  = length(vX);  % number of grid points along x-axis
    nY  = length(vY);
    
    xmin   = min(vX); % xmin
    fDx    = max(vX) - xmin; % xmax - xmin
    ymin   = min(vY); % ymin
    fDy    = max(vY) - ymin; % ymax - ymin
    xlimit = [xmin xmin+fDx];
    ylimit = [ymin ymin+fDy];
    xbox   = xlimit([1 1 2 2 1]);
    ybox   = ylimit([1 2 2 1 1]);
    vFlist = zeros(nF,4);  % the fault database [x1,y1,x2,y2]
    vFpoint = [];         % the points placed on each fault equidistantly
    vHj = [];              % the witdth for the Gaussian kernel for each fault
    % generate nF faults inside a rectangular box
    for n = 1:nF
        r = d*sqrt(rand()^(1/(1-q)) - 1); % the length of a fault drawn from a distribution with parameters d and q
        phi = 2*pi*rand();                % the uniform random orientation of faults
        % random position in a rectangular box
        x1 = xmin + fDx*rand();           % x-coordinate of the starting point of a fault
        y1 = ymin + fDy*rand();           % y-coordiante
        x2 = x1 + r*cos(phi);             % the end point of a fault
        y2 = y1 + r*sin(phi);
        % if the end point of the fault is outside the rectangle then find
        % a point on a boundary and make it the end point
        [xi,yi] = polyxpoly([x1 x2],[y1 y2],xbox,ybox); % check the interstion of a fault line with the boundary
        if ~isempty(xi) || ~isempty(yi)    % if there is an intersection replace the end (x2,y2) with the itersection point on the boundary
            x2 = xi;
            y2 = yi;
        end
        
        vFlist(n,:) = [x1, y1, x2, y2];   % add the fault into the list
        %dr = 1/sqrt((x1-x2)^2+(y1-y2)^2);
        dr = d/sqrt((x1-x2)^2+(y1-y2)^2);
        %disp(dr)
        vxn = (x1:dr*cos(phi):x2)';       % add equidistant points along the fault
        vyn = (y1:dr*sin(phi):y2)';
        vFpoint = [vFpoint; [vxn, vyn]]; % add these point to the total list
        vHj = [vHj; 3.5*ones(size(vxn))]; % for each fault compute the bandwidth for the Gaussian kernel
    end
    %disp(vHj)

    fSigma = 1.0;
    lsize  = max(nX,nY);
    nMax   = nextpow2(lsize);
    mUfb   = fbsurf(nMax,fSigma,fH,true);
    fUmax  = max(max(mUfb(1:nY,1:nX)));
    fUmin  = min(min(mUfb(1:nY,1:nX)));
    mUfb   = (mUfb(1:nY,1:nX) - fUmin)/(fUmax - fUmin);   % scale it between 0 and 1
    %mUfb   = mUfb(1:nY,1:nX) - fUmin;

    mUfault = zeros(nY,nX);
    nJ = length(vFpoint(:,1));
    for iy = 1:nY
        for ix = 1:nX
            for j = 1:nJ
                r = sqrt((vX(ix) - vFpoint(j,1))^2 + (vY(iy) - vFpoint(j,2))^2);
                fZj = gauss_kernel(r,vHj(j));
                mUfault(iy,ix) = mUfault(iy,ix) + fZj;
            end
        end
    end
    fUmax = max(max(mUfault));
    fUmin = min(min(mUfault));
    mUfault = (mUfault - fUmin)/(fUmax - fUmin);   % scale it between 0 and 1
    
    mUtmp = fP*mUfb + (1-fP)*mUfault;
    
    fUmax = max(max(mUtmp));
    fUmin = min(min(mUtmp));
    mUtmp = (mUtmp - fUmin)/(fUmax - fUmin);   % scale it between 0 and 1
    
    indx  = mUtmp > fI0;                                     % consider only top (1-fL)*100% to create regions with 0 rate
    mU    = zeros(nY,nX);
    mU(indx) = mUtmp(indx);
    if bNormFactor
        mU = fNormFactor.*mU;
    end
end
