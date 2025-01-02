function rate_cum = etas2d8p_cum(fT,vPar,vCat,ETAS,Bckgrnd_est)
%
%   
%   using the normalization by Zhuang et al (2005)
%
%   vMu    - the background rate
%   vPar   - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%   vCat   - earthquake times and magnitudes with the subtracted reference magnitude m_0 (vTM(:,2) - m_0)
%   fTs    - the start time for the target window
%   fTe    - the end time for the target window
%   nJs    - the first index of the event in the target window
%   nJe    - the last index of the event in the target window
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 12 December 2021
%   ...
%   version: 1.0.0, 12 December 2021
%
    nJs     = ETAS.nJs;
    nJ      = find(vCat(:,1) < fT,1,'last'); % find the last index for which vCat(:,1) < vT(jt)
    vexpalf = exp(vPar(3).*vCat(1:nJ,4)); % exp(alpha*m_i) assumes that m_0 is subtracted
    vexpgam = exp(vPar(8).*vCat(1:nJ,4)); % exp(gamma*m_i)
    vtc     = vCat(1:nJ,1)./vPar(4);
    fTsc    = ETAS.fTs/vPar(4);            % fTs/c
    fTec    = fT/vPar(4);                  % fTe/c
    p1      = 1.0 - vPar(5);               % 1 - p
    nJs1    = nJs - 1;
    
    spat     = etas2d8p_spat_prod(ETAS.vRk,ETAS.vRk2,ETAS.vDk,nJ,ETAS.inR,vPar(6),vPar(7),vexpgam);
    rate_cum = sum(vexpalf(1:nJs1).*spat(1:nJs1).*((fTsc - vtc(1:nJs1) + 1).^p1 - (fTec - vtc(1:nJs1) + 1).^p1)); % all events in the catalog for t < Ts
    rate_cum = vPar(2).*(rate_cum + sum(vexpalf(nJs:nJ).*spat(nJs:nJ).*(1 - (fTec - vtc(nJs:nJ) + 1).^p1)));  %
    
    vBgProb_est = Bckgrnd_est.vBgProb;
    fUquadR  = quad_background(ETAS.inR,nJ,vBgProb_est,ETAS.vHj,ETAS.vRk,ETAS.vRk2,ETAS.vDk,ETAS.nNk)/(fT - ETAS.fT0);
    rate_cum = vPar(1)*fUquadR*(fT - ETAS.fTs) + rate_cum;
end

