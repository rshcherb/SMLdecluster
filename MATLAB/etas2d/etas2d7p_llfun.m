function fLle = etas2d7p_llfun(vPar,vCat,ETAS)
%
%   The log-likelihood function for the spatial ETAS conditional rate using the normalization by Zhuang et al (2005)
%
%   vPar   - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7)
%   vCat   - earthquake times and magnitudes with the subtracted reference magnitude m_0 (vTM(:,2) - m_0)
%   ETAS   - structure with fields
%       vU     - the interpolated background rate u(x,y) at the epicentres
%       fUquadR - the spatial integral of the rate function u(x,y) over the region R
%       fTs    - the start time for the target window [T_s,T_e]
%       fTe    - the end time for the target window [T_s,T_e]
%       nJs    - the first index of the event in the target window [T_s,T_e]
%       nJe    - the last index of the event in the target window [T_s,T_e]
%       inR    - logical indices of events in the target region R during [T0, Te]
%       linTR  - linear indices of each event in the target region R and during [T_s,T_e]
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 5 May 2021
%   ...
%   version: 1.0.0, 31 May 2021
%
    vexpalf = exp(vPar(3).*vCat(1:ETAS.nJe,4)); % assumes that m_0 is subtracted
    vtc    = vCat(1:ETAS.nJe,1)./vPar(4);
    d2     = vPar(6)^2;
    vx     = vCat(1:ETAS.nJe,3);                % x/longitude
    vy     = vCat(1:ETAS.nJe,2);                % y/latitude
    fFac   = vPar(2)*(vPar(5)-1)/vPar(4)*(vPar(7)-1)/(pi*d2); % A*(p-1)/c*(q-1)/(\pi*d^2)
    sumlog_gif = 0.0;
    for j = (ETAS.linTR')         % the earthquakes during the target time interval t_j \in [Ts, Te] and inside the target region (x_j,y_j) \in R
        s = vPar(1)*ETAS.vU(j) + fFac*sum( vexpalf(1:(j-1)).*((1 + vtc(j) - vtc(1:(j-1))).^(-vPar(5)) ).*...
            ((1.0 + ( (vx(j) - vx(1:(j-1))).^2 + (vy(j) - vy(1:(j-1))).^2 )/d2).^(-vPar(7))) ); % \mu + sum
        if s > 1e-25
            sumlog_gif = sumlog_gif + log(s);
        else
            sumlog_gif = sumlog_gif - 50;
        end
    end
    
    fLle = -etas2d7p_prod(vPar,vCat,ETAS) + sumlog_gif; %
    %disp(vPar)
    %disp(fLle)
end
