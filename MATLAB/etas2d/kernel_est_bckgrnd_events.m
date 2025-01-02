function [vUest, vUtot, vClust] = kernel_est_bckgrnd_events(fT0,fT,vCat,vHj,vBgProb)
%
%   Computes the background rate using the variabe kernel estimate (Zhuang, EPS, 2011, p.208, Eq. (10); Zhuang et al., JGR 2005, p.4, Eq. (20))
%   This is done at each location of the events in the catalogue
%
%   fT0    - time start of the catalogue
%   fT     - time end of the study interval
%   vCat   - event catalogue
%   vHj    - the bandwidth (standard deviation) from spat_bandwidth.m to be used in the Gaussian kernel
%   vBgProb - the probability that event j is a background/spontaneous event
%
%   vUest  - variable kernel estimate of the spatial background rate at the earthquake epicentres
%   vUtot  - total seismicity rate
%   vClust - clustering coefficient
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 19 May 2021
%   ...
%   version: 1.0.0, 14 June 2021
%
    nJ    = find(vCat(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= vT
    vUtot = zeros(nJ,1);
    vUest = zeros(nJ,1);
    for i = 1:nJ    % for the location of each event i compute the contribution to the background rate from all events
        for j = 1:nJ
            %fZj = mvnpdf([vCat(i,3) - vCat(j,3),vCat(i,2) - vCat(j,2)],[0 0],[vHj(j) 0; 0 vHj(j)]);
            r = sqrt((vCat(i,3) - vCat(j,3))^2 + (vCat(i,2) - vCat(j,2))^2);
            fZj = gauss_kernel(r,vHj(j));
            vUtot(i)  = vUtot(i) + fZj;
            vUest(i) = vUest(i) + vBgProb(j)*fZj;
        end
    end
    vClust = 1 - vUest./vUtot;
    vUtot  = vUtot./(fT - fT0);
    vUest  = vUest./(fT - fT0);
end
