function vRate = gif_etas2d8p(vT,fX,fY,fU,vPar,d2,fFac,vCat)
%
%   Conditional ground intensity rate for the spatial ETAS model with 8 parameters
%
%   vT   - times at which to compute the ground intensity function
%   fX   - X coordinate (longitude)
%   fY   - Y coordinate (latitude)
%   fU   - the background rate u(x,y) at the point (fX, fY): 
%          the point (fX,fY) can be either the location of events from the vCat, or any arbitrary location within the region
%   vPar - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%   d2   - d^2 (vPar(6)^2)
%   fFac - A*(p-1)/c*(q-1)/(\pi*d^2)
%   vCat - earthquake catalogue with the subtracted reference magnitude m_0 (vTM(:,2) - m_0)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 16 May 2021
%   ...
%   version: 1.0.0, 14 June 2021
%
    nT    = length(vT);
    vRate = zeros(1,nT);
    for jt = 1:nT
        nJ = find(vCat(:,1) < vT(jt),1,'last'); % find the last index for which vCat(:,1) < vT(jt)
        vRate(jt) = vPar(1)*fU + ...
                  fFac*sum( exp((vPar(3)-vPar(8)).*vCat(1:nJ,4)).*...
                          (((vT(jt) - vCat(1:nJ,1))./vPar(4) + 1).^(-vPar(5))).*...
                          ((1 + ( (fX - vCat(1:nJ,3)).^2 + (fY - vCat(1:nJ,2)).^2 )./(d2.*exp(vPar(8).*vCat(1:nJ,4)))).^(-vPar(7))) ); % \mu + sum
    end
end

