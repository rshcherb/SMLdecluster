function [vX, vMagHist, vMagCum] = gr_fit_plot(vMag,vPar,vParErr,fMc,fDm,varargin)
%
%   Function to plot the frequency-magnitude statistics, GR scaling, Bath's law
%
%   vMag    - list a magnitudes to bin and construct the distribution
%   vPar    - the parameters of GR scaling
%   vParErr - the errors of the parameters
%   fMc     - the lower magnitude cutoff used to estimate parameters
%   fDm     - magnitude binning
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 20 November 2019
%   ...
%   version: 1.0.0, 20 November 2019
%
    fAlpha   = 0.05;
    bPoisson = false;
    bBath    = false;
    bPlotAll = true; bPlotHist = false; bPlotCumul = false;
    fXmin    = floor(min(vMag)/fDm)*fDm;
    fXmax    = ceil(max(vMag)/fDm)*fDm; % fix the upper limit
    for k = 1:length(varargin)
        if strcmp('Xmin',varargin{k})
            fXmin  = varargin{k+1};
        end
        if strcmp('Xmax',varargin{k})
            fXmax  = varargin{k+1};
        end
        if strcmp('Histogram',varargin{k})
            bPlotHist = true;
            bPlotAll  = false;
        end
        if strcmp('Cumulative',varargin{k})
            bPlotCumul = true;
            bPlotAll   = false;
        end
        if strcmp('PoissonBounds',varargin{k})
            bPoisson = true;
        end
        if strcmp('BathLaw',varargin{k})
            bBath = true;
        end
        if strcmp('Alpha',varargin{k})
            fAlpha = varargin{k+1};
        end
    end
    
    vX = fXmin:fDm:fXmax;   % centres of the bins, total N, between Xmin and Xmax 
    %                                                            vX(1)        vX(k)
    vXedges = (fXmin-0.5*fDm):fDm:(fXmax+0.5*fDm); % bin edges: |_____|_..._|_____|_..._|_____|, total N+1
    %                                                           1     2     k    k+1
    % compute the histogram
    vMagHist = histcounts(vMag,vXedges);           % histogram of the earthquake magnitudes
    nBin = length(vMagHist);                       % number of bins, N

    % complimentary cumulative distribution
    vMagCum = zeros(1,nBin);
    fSum = 0.0;
    for j = nBin:-1:1
        fSum = fSum + vMagHist(j);
        vMagCum(j) = fSum;
    end
    indx = vMagHist > 0; % use only nonzero bins
    vX = vX(indx);
    vMagHist = vMagHist(indx);
    vMagCum = vMagCum(indx);
    
    if bPlotAll || bPlotCumul
        semilogy(vX,vMagCum,'s','MarkerEdgeColor','#D95319');
        hold on
        semilogy(vX,10.^(-vPar(1).*vX + vPar(2)),'Color','black');
    end
    if bPlotAll || bPlotHist
        hold on
        semilogy(vX,vMagHist,'s','MarkerEdgeColor','#4DBEEE','MarkerFaceColor','#4DBEEE','MarkerSize',4);
        %bar(vX,vMagHist,0.95,'FaceColor','#4DBEEE','EdgeColor','#4DBEEE');
        if bPoisson
            % compute Poisson confidence interval for the histogram
            [vXpdf, vPdfGR, vPdfErrLo, vPdfErrHi] = gr_poisson(vX,vMagHist,fMc,fAlpha);
            semilogy(vXpdf,vPdfGR,'-','Color','#7E2F8E');
            semilogy(vXpdf,vPdfErrLo,'--','Color','#7E2F8E');
            semilogy(vXpdf,vPdfErrHi,'--','Color','#7E2F8E');
        end
    end
    set(gca,'XMinorTick','on');
    set(gca,'Layer','top');
    sposx = 0.7;
    sposy = 0.9;
    dy = 0.05;
    str = sprintf('b = %.3f \xB1 %.3f',vPar(1),vParErr(1));
    text(sposx,sposy,str,'Units','normal');
    str = sprintf('a = %.3f \xB1 %.3f',vPar(2),vParErr(2));
    sposy = sposy-dy;
    text(sposx,sposy,str,'Units','normal');
    xlabel('magnitude, m');
    ylabel('N(\geq m)');
    if bBath
        % compute the inferred largest earthquake (modified Bath law)
        fMstar = gr_bath(fMc,vPar,vParErr);
        % plot m* magnitude
        plot(fMstar,1,'ks','MarkerSize',10,'Linewidth',1.5);
        plot(fMstar,1,'k+','MarkerSize',8,'Linewidth',1.5);
        sposy = sposy-dy;
        text(sposx,sposy,sprintf('m* = %.3f',fMstar),'Units','normal');
    end
    % plot mc cutoff line
    plot([fMc, fMc],[0.8 max(vMagCum)],'--');
    axis tight;
    xlim([fXmin, Inf]);
    ylim([0.8, Inf]);
end
