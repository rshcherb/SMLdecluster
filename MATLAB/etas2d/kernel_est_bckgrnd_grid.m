function [mUest, mUtot, mClust] = kernel_est_bckgrnd_grid(fT0,fT,vX,vY,vCat,vHj,vBgProb)
%
%   Computes the background rate using the variabe kernel estimate (Zhuang, EPS, 2011, p.208, Eq. (10); Zhuang et al., JGR 2005, p.4, Eq. (20))
%   This is done for each point of the rectangular grid (vX,vY)
%
%   fT0    - time start of the catalogue
%   fT     - time end of the study interval
%   vX     - vector for X coordinates of the renctanglular grid
%   vY     - vector for Y coordinates of the renctanglular grid
%   vCat   - event catalogue
%   vHj    - the bandwidth (standard deviation) from spat_bandwidth.m to be used in the Gaussian kernel
%   vBprob - the probability that event j is a background/spontaneous event
%
%   mUest  - variable kernel estimate of the spatial background rate
%   mUtot  - total seismicity rate
%   mClust - clustering coefficient
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 19 May 2021
%   ...
%   version: 1.0.0, 14 June 2021
%
    nX    = length(vX);
    nY    = length(vY);
    mUtot = zeros(nY,nX);
    mUest = zeros(nY,nX);
    nJ    = find(vCat(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= vT(j)
    for iy = 1:nY
        for ix = 1:nX
            for j = 1:nJ
                %fZj = mvnpdf([vX(ix) - vCat(j,3),vY(iy) - vCat(j,2)],[0 0],[vHj(j) 0; 0 vHj(j)]);
                r = sqrt((vX(ix) - vCat(j,3))^2 + (vY(iy) - vCat(j,2))^2);
                fZj = gauss_kernel(r,vHj(j));
                mUtot(iy,ix)  = mUtot(iy,ix) + fZj;
                mUest(iy,ix) = mUest(iy,ix) + vBgProb(j)*fZj;
            end
        end
    end
    mClust = 1 - mUest./mUtot;
    mUtot  = mUtot./(fT - fT0);
    mUest  = mUest./(fT - fT0);
end
