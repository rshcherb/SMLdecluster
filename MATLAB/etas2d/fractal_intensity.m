function mU = fractal_intensity(nX,nY,fH,varargin)
%
%   Generates a background earthquake rate u(x,y) on a rectangular grid based on the 2D fractal landscape with Hurst exponent H
%
%   nX, nY - the size of the grid (256,256)
%   fH     - the Hurst exponent: (0, 1)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 17 November 2014
%   ...
%   version: 2.0.0, 18 May 2021
%
    bNormFactor = false;
    for k = 1:length(varargin)
        if strcmp('NormFactor',varargin{k})
            bNormFactor = true;
            fNormFactor = varargin{k+1};
        end
    end
    
    fSigma = 1.0;
    lsize  = max(nX,nY);
    nMax   = nextpow2(lsize);
    tmpU   = fbsurf(nMax,fSigma,fH,true);
    tmpMax = max(max(tmpU(1:nY,1:nX)));
    tmpMin = min(min(tmpU(1:nY,1:nX)));
    %disp([tmpMin, tmpMax])
%
%     % normalize to be between [0, 1]
%     mU = (tmpU(1:nY,1:nX) - tmpMin)/(tmpMax - tmpMin); % normalize to be between 0 and 1
%
    % 
    tmpU   = (tmpU(1:nY,1:nX) - tmpMin)/(tmpMax - tmpMin);   % scale it between 0 and 1
    indx   = tmpU > 0.3;                                     % consider only top 70% to create regions with 0 rate
    mU     = zeros(nY,nX);
    mU(indx) = tmpU(indx);
%
%     % clip below zero
%     tmpU = tmpU(1:nY,1:nX);
%     indx = tmpU > 0;  % consider only cells above 0
%     mU = zeros(nY,nX);
%     mU(indx) = tmpU(indx);
%
    if bNormFactor
        mU = fNormFactor.*mU;
    end
end
