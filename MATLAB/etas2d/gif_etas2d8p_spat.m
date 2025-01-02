function fRate = gif_etas2d8p_spat(fTs,fTe,nJs,nJe,fX,fY,fMuU,vPar,vCat)
%
%   Conditional ground intensity rate for the spatial ETAS model with 8
%   parameters integrated from fTs to fTe
%
%   fTs, fTe - 
%   nJs, nJe - 
%   fX   - X coordinate (longitude)
%   fY   - Y coordinate (latitude)
%   fMuU - the background rate \mu*u(x,y) at the point (fX, fY): 
%          the point (fX,fY) can be either the location of events from the vCat, or any arbitrary location within the region
%   vPar - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%   vCat - earthquake catalogue with the subtracted reference magnitude m_0 (vTM(:,2) - m_0)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 24 December 2021
%   ...
%   version: 1.0.0, 24 December 2021
%
    vexpalf = exp(vPar(3).*vCat(1:nJe,4)); % exp(alpha*m_i) assumes that m_0 is subtracted
    %vexpgam = exp(vPar(8).*vCat(1:nJe,4)); % exp(gamma*m_i)
    vtc     = vCat(1:nJe,1)./vPar(4);
    fTsc    = fTs/vPar(4);            % fTs/c
    fTec    = fTe/vPar(4);            % fTe/c
    p1      = 1.0 - vPar(5);          % 1 - p
    nJs1    = nJs - 1;
    d2      = vPar(6)^2;              % d^2
    fFac    = vPar(2)*(-p1)/(vPar(4))*(vPar(7)-1)/(pi*d2); %   fFac = A*(p-1)/c*(q-1)/(\pi*d^2)
    
    spat = ((1 + ( (fX - vCat(1:nJe,3)).^2 + (fY - vCat(1:nJe,2)).^2 )./(d2.*exp(vPar(8).*vCat(1:nJe,4)))).^(-vPar(7)));
    lamb_sum = sum(vexpalf(1:nJs1).*spat(1:nJs1).*((fTsc - vtc(1:nJs1) + 1).^p1 - (fTec - vtc(1:nJs1) + 1).^p1)); % all events in the catalog for t < Ts
    lamb_sum = fFac.*(lamb_sum + sum(vexpalf(nJs:nJe).*spat(nJs:nJe).*(1 - (fTec - vtc(nJs:nJe) + 1).^p1))); %   

    fRate = fMuU*(fTe - fTs) + lamb_sum;
end

