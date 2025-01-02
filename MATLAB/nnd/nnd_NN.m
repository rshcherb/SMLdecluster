function NND = nnd_NN(vCat,Par,options)
%
%   Performs NND computation and classification of events without any
%   thresholds on \eta
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.1.0, 3 December 2024
%
    arguments
        vCat double
        Par struct
        options.NN double = 1
        options.FMD char = 'GR'            % frequency-magnitude distribution to rescale \eta, T, and R: 'GR', 
        options.Distance char = 'Radial'   % distance between events: 'Radial', 'Euclidean', 'Euclidean3D'   
    end

    %cp = 6.0*86400.0; %*365; % p-wave velocity in km/s
    vMag = vCat(:,4);
    if strcmp (options.FMD,'GR')
        b = Par.b; df = Par.df; q = Par.q; rmin = Par.r_min; rmax = Par.r_max; p = 1 - q;
        FMDi = 10.^(-b*vMag);
        FMDi_T = 10.^(-q*b*vMag);
        FMDi_R = 10.^(-p*b*vMag);
    elseif strcmp(options.FMD,'Exp')
        beta = Par.beta; m_min = Par.m_min; df = Par.df; q = Par.q; rmin = Par.r_min; rmax = Par.r_max; p = 1 - q;
        FMDi = exppdf((vMag - m_min),1.0/beta);
        FMDi_T = FMDi.^q;
        FMDi_R = FMDi.^p;
    elseif strcmp(options.FMD,'TapTapPareto')
        vMom = mag2moment(vMag);
        M_min = Par.M_min; alpha1 = Par.alpha1; M_cm1 = Par.M_cm1; alpha2 = Par.alpha2; M_cm2 = Par.M_cm2; w = Par.w;
        df = Par.df; q = Par.q; rmin = Par.r_min; rmax = Par.r_max; p = 1 - q;
        FMDi_mom = taptap_pareto_pdf(vMom,M_min,alpha1,M_cm1,alpha2,M_cm2,w);
        FMDi = pdf_moment2mag(vMag,FMDi_mom); % convert the pdf from moment domain to magnitude domain
        FMDi_T = FMDi.^q;
        FMDi_R = FMDi.^p;
    elseif strcmp(options.FMD,'TapParetoPareto')
        vMom = mag2moment(vMag);
        M_min = Par.M_min; alpha = Par.alpha; M_cm = Par.M_cm; beta = Par.beta; w = Par.w;
        df = Par.df; q = Par.q; rmin = Par.r_min; rmax = Par.r_max; p = 1 - q;
        FMDi_mom = tap_pareto_pareto_pdf(vMom,M_min,alpha,M_cm,beta,w);
        FMDi = pdf_moment2mag(vMag,FMDi_mom); % convert the pdf from moment domain to magnitude domain
        FMDi_T = FMDi.^q;
        FMDi_R = FMDi.^p;
    end

    if strcmp(options.Distance,'Radial')
        distance_ij = @(lati,loni,latj,lonj) radialdist(lati,loni,latj,lonj);
    elseif strcmp(options.Distance,'Euclidean')
        distance_ij = @(yi,xi,yj,xj) sqrt((xi - xj)^2 + (yi - yj)^2);
    end

    nEq     = size(vCat,1);
    mEta_ij = Inf(nEq,nEq);
    mR_ij   = Inf(nEq,nEq);
    mT_ij   = Inf(nEq,nEq);
    for i = 1:nEq-1
        for j = (i+1):nEq
            t_ij = vCat(j,1) - vCat(i,1); % j is a subsequent event
            r_ij = distance_ij(vCat(i,2),vCat(i,3),vCat(j,2),vCat(j,3));
            %if t_ij > 0.0 && r_ij > rmin && r_ij <= cp*t_ij % 
            %if t_ij > 0.0 && r_ij >= rmin && r_ij <= rmax % excludes events that are either occurred at the same time
            if t_ij > 0.0 && r_ij > rmin && r_ij <= rmax % excludes events that are either occurred at the same place or time
                mEta_ij(i,j) = t_ij*(r_ij)^df*FMDi(i);
                mT_ij(i,j)   = t_ij*FMDi_T(i);
                mR_ij(i,j)   = (r_ij)^df*FMDi_R(i);
            end
        end
    end
    NND.vEta    = Inf(nEq,options.NN);   % nearest-neighbour proximity
    NND.vT      = Inf(nEq,options.NN);   % rescaled time
    NND.vR      = Inf(nEq,options.NN);   % rescaled distance
    NND.vP_indx = ones(nEq,1);   % indeces of parent events in vCat
    NND.vP_mag  = zeros(nEq,1);  % magnitudes of parent events
    NND.vCh_mag = zeros(nEq,1);  % magnitudes of child events
    NND.vDmag   = zeros(nEq,1);  % magnitude differnce between the parent and child events
    NND.vDt     = Inf(nEq,1);    % the time difference between the parent and child events
    NND.vDist   = Inf(nEq,1);    % the distance in km between the parent and child events
    NND.vClust  = zeros(nEq,1);
    NND.vP_mag(1)  = vCat(1,4);  % the first event is a singleton and does not have a parent
    NND.vCh_mag(1) = vCat(1,4);
    for j = 2:nEq
        %[vEta, pind] = min(mEta_ij(1:j-1,j));
        %[vEta, indx] = sort(mEta_ij(1:j-1,j));
        [vEta, indx] = mink(mEta_ij(1:j-1,j),options.NN);
        vEta = vEta';
        pind = indx(1);
        n = length(vEta);
        if n >= options.NN
            NND.vEta(j,:)    = vEta(1:options.NN);
            NND.vT(j,:)      = mT_ij(indx(1:options.NN),j);
            NND.vR(j,:)      = mR_ij(indx(1:options.NN),j);
        else
            NND.vEta(j,1:n)    = vEta(1:n);
            NND.vT(j,1:n)      = mT_ij(indx(1:n),j);
            NND.vR(j,1:n)      = mR_ij(indx(1:n),j);
        end
        NND.vDt(j)     = vCat(j,1) - vCat(pind,1);
        NND.vDist(j)   = distance_ij(vCat(pind,2),vCat(pind,3),vCat(j,2),vCat(j,3));
        NND.vP_indx(j) = pind;
        NND.vP_mag(j)  = vCat(pind,4);
        NND.vCh_mag(j) = vCat(j,4);
        NND.vDmag(j)   = NND.vP_mag(j) - NND.vCh_mag(j);
        %NND.vClust(j)  = 1;
    end
end

