function etas_prod = etas2d7p_prod(vPar,vCat,ETAS)
%
%   The log-likelihood function for the spatial ETAS conditional rate
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
%   version: 1.0.0, 5 May 2021
%   ...
%   version: 1.0.0, 31 May 2021
%
    nJs     = ETAS.nJs;
    nJe     = ETAS.nJe;
    vexpalf = exp(vPar(3).*vCat(1:nJe,4)); % exp(alpha*m_i) assumes that m_0 is subtracted
    vtc     = vCat(1:nJe,1)./vPar(4);
    fTsc    = ETAS.fTs/vPar(4);            % fTs/c
    fTec    = ETAS.fTe/vPar(4);            % fTe/c
    p1      = 1.0 - vPar(5);               % 1 - p
    nJs1    = nJs - 1;
    
    spat = etas2d7p_spat_prod(ETAS.vRk,ETAS.vRk2,ETAS.vDk,nJe,ETAS.inR,vPar(6),vPar(7));
    etas_prod = sum(vexpalf(1:nJs1).*spat(1:nJs1).*((fTsc - vtc(1:nJs1) + 1).^p1 - (fTec - vtc(1:nJs1) + 1).^p1)); % all events in the catalog for t < Ts
    etas_prod = vPar(2).*(etas_prod + sum(vexpalf(nJs:nJe).*spat(nJs:nJe).*(1 - (fTec - vtc(nJs:nJe) + 1).^p1))); %   
    
    etas_prod = vPar(1)*ETAS.fUquadR*(ETAS.fTe-ETAS.fTs) + etas_prod;
end

