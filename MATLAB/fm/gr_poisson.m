function [vXpdf, vPdfGR, vPdfErrLo, vPdfErrHi] = gr_poisson(vXbin,vBinData,fMc,fAlpha)
%
%   Computes the poisson confidence intervals for the histogram
%
%   vXbin    - coordinates of x bins
%   vBinData - histogram
%   fMc      - lower magnitude cutoff
%   fAlpha   - confidence level: 0.05 (95%)
%
% Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
% version: 1.0.0, October 20, 2014
%
    xindx = vXbin >= fMc;
    vXpdf = vXbin(xindx);
    [histmodel, dev, stats] = glmfit(vXpdf,vBinData(xindx),'poisson','link','log');
    [vPdfGR, vPdfLo, vPdfHi] = glmval(histmodel,vXpdf,'log',stats,'confidence',1.0-fAlpha);

    % compute Poisson confidence intervals
    vDataErrLow = poissinv(0.5*fAlpha,vPdfGR);
    vDataErrUp  = poissinv(1.0-0.5*fAlpha,vPdfGR);

    vPdfErrLo = vDataErrLow;
    vPdfErrHi = vDataErrUp;
end
