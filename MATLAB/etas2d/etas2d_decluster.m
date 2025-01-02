function [vBgDeclstr, indxBg] = etas2d_decluster(fT,vCat,vBgProb)
%
%   Outputs background/spontaneous events and the index logical array:  (Zhuang et al., JGR 2005, p.4, Eq. (18); Zhuang, EPS 63 (2011), p.208)
%
%   fT         - the maximu time for the catalogue 
%   vCat       - event catalogue
%   vBgProb    - the background probability for each event
%   vBgDeclstr - the declustered events
%   indxBg     - the logical index: 1 - background events, 0 - triggered events (aftershocks)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 1 June 2021
%   ...
%   version: 1.0.0, 1 June 2021
%
    nJ = find(vCat(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= vT(j)
    %disp([length(vCat(:,1)), nJ])
    indxBg = false(nJ,1);
    n = 0;
    for i = 1:nJ
        if  rand() < vBgProb(i)
            n = n + 1;
            vBgDeclstr(n,:) = vCat(i,:);  % background events from stochastic declustering
            indxBg(i) = true;
        end
    end
end
