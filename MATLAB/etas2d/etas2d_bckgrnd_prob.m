function vBgProb = etas2d_bckgrnd_prob(vCat,vU,vPar,fT,varargin)
%
%   Computes the probability that event i is a background/spontaneous event:  (Zhuang et al., JGR 2005, p.4, Eq. (18); Zhuang, EPS 63 (2011), p.208)
%
%   vCat - event catalogue: magnitudes must be m - m_0 
%   vU   - the value of the background rate u(x,y) at the location of each event i up to fT
%   vPar - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%   fT   - the end of the study time interval
%   vBgProb - the background probability for for the locations of each event
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 18 May 2021
%   ...
%   version: 1.1.0, 29 December 2021
%
    sPointProc = 'ETAS2D8P'; % 
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
    end
    
    d2     = vPar(6)^2;
    fFac   = vPar(2)*(vPar(5)-1)/vPar(4)*(vPar(7)-1)/(pi*d2); % A*(p-1)/c*(q-1)/(\pi*d^2)
    nJ     = find(vCat(:,1) <= fT,1,'last');                  % find the last index for which vCat(:,1) <= fT   
    vBgProb = zeros(nJ,1);
    if strcmp(sPointProc,'ETAS2D7P')
        for i = 1:nJ
            vBgProb(i) = vPar(1)*vU(i)./gif_etas2d7p(vCat(i,1),vCat(i,3),vCat(i,2),vU(i),vPar,d2,fFac,vCat);
        end
    elseif strcmp(sPointProc,'ETAS2D8P')
        for i = 1:nJ
            vBgProb(i) = vPar(1)*vU(i)./gif_etas2d8p(vCat(i,1),vCat(i,3),vCat(i,2),vU(i),vPar,d2,fFac,vCat);
        end
    end
    %disp([min(vBprob),mean(vBprob),median(vBprob),max(vBprob)])
end
