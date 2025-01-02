function res = etas2d7p_spat_prod(vRk,vRk2,vDk,nJe,inR,d,q)
%
%   Computes the spatial integral for the productivity part of he likelihood function by considering the
%   subdivision of the region R into radial segments with each event being a centre
%
%   vRk  - the mean radius of each radial segment: size [nEQ, nVert_new]
%   vDk  - the angle of each radial segement in radians
%   inR  - the logical indeces of events in the catalogue that are in R
%   d, q - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 5 May 2021
%   ...
%   version: 1.0.0, 31 May 2021
%
    res   = zeros(nJe,1);
    f12pi = 1.0/(2*pi);
    f1d2  = 1.0/(d^2);
    q1    = 1.0 - q;
    for n = 1:nJe
        if inR(n)      % for events inside the target region R
            res(n) = f12pi*sum(vDk(n,:).*(1.0 - (1.0 + vRk(n,:).^2*f1d2).^q1));
        else
            resR = f12pi*sum(vDk(n,:).*(1.0 - (1.0 + vRk(n,:).^2*f1d2).^q1));
            resR2 = f12pi*sum(vDk(n,:).*(1.0 - (1.0 + vRk2(n,:).^2*f1d2).^q1));
            res(n) = resR2 - resR;
        end
    end
end

