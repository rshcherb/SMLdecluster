function etas2d_sequence(Parent,beta1,vPar,fMmin,fM0,fTmax,nNmax,sPointProc,vBndr)
%
%   Parent - structure of paraent event with fields: time, lat, lon, mag, rank
%   beta1  - 1/fBeta: \beta = log(10)*b, b-value for the GR scaling law
%   vPar   - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%   fMmin  - the minimum magnitude cutoff it is assumed that fMmin = fMc - fDm
%   fM0    - the reference magnitude
%   fTmax  - the maximum time interval to use for aftershocks
%   nNmax  - the absolute mixumum number of events to simulate
%   vEqCat - global variable for the event catalogue
%   vEqNum - global variable for the total event count
%
%        vTM            vTMout
%   fT0 ------- fTs ---------------- fTe
%                future time window
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 April 2021
%   ...
%   version: 1.0.0, 10 May 2021
%
    global vEqCat;     % earthquake catalogue
    global nEqNum;     % earthquake counter

    nProd  = poissrnd(vPar(2)*exp(vPar(3)*(Parent.mag - fM0)));
    asrank = Parent.rank + 1;
    d = vPar(6);
    if strcmp(sPointProc,'ETAS2D8P')
        d = vPar(6)*exp(0.5*vPar(8)*(Parent.mag - fM0));     % for Model 5
    end
    for n = 1:nProd
        fAsTime = Parent.time + mol_inverse(vPar(4),vPar(5)); % for each aftershock the time is drawn from the normalized distribution 
        if fAsTime <= fTmax && nEqNum <= nNmax
            fMag  = fMmin + exprnd(beta1);          % exponential magnitude
            [lon, lat] = dr_rad_inverse(Parent.lon,Parent.lat,d,vPar(7)); % (Px,Py,d,q) power-law distributed random position from the parent event epicentre
            if lon > vBndr(1) && lon < vBndr(2) && lat > vBndr(3) && lat < vBndr(4) % ony accept events in vBound box
                nEqNum = nEqNum + 1;
                vEqCat(nEqNum,:) = [fAsTime,lat,lon,fMag,asrank]; % update the catalogue
                if fMag > fM0
                    Parent_new.time = fAsTime;
                    Parent_new.lon  = lon;
                    Parent_new.lat  = lat;
                    Parent_new.mag  = fMag;
                    Parent_new.rank = asrank;
                    etas2d_sequence(Parent_new,beta1,vPar,fMmin,fM0,fTmax,nNmax,sPointProc,vBndr);
                end
            end
        end
    end
end
