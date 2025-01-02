function NND = nnd_decluster_zbz2020(vCat,NND,NNDpar,vEtaThresh,sDistance)
%
%   Performs NND declustering based on Zaliapin and Ben-Zion, JGR (2020) doi: 10.1029/2018JB017120
%   
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 12 October 2024
%   ...
%   version 1.2.0, 2 December 2024
%
    nEq = size(vCat,1);
    NND = nnd_decluster_eta(NND,vEtaThresh,NNDpar.dt_max,NNDpar.dr_max); % decluster the catalog based on the original NND algorithm

    indx = NND.vClust == 0; % identify background events based on the original NND algorithm with vEtaThresh
    C0 = vCat(indx,:);      % the catalog of only background events
    N0 = length(C0(:,1));
    t0 = min(vCat(:,1));    % the minumum time oin the original vCat catalog
    t1 = max(vCat(:,1));
    M = NNDpar.M;
    A0 = 10^(NNDpar.alpha0);
    kappa = zeros(nEq,M);
    for k = 1:M
        Ck = C0;
        Ck(:,4) = Ck(randperm(N0)',4);             % random permute magnitudes
        Ck(:,1) = sort(t0 + (t1 - t0)*rand(N0,1)); % replace with sorted random times
        for j = 1:nEq
            kappa(j,k) = nn_prox(vCat(j,:),Ck,NNDpar,sDistance); % compute the NN proximity for each event in the original catalog
        end
    end
    log10alpha = zeros(1,nEq);
    for j = 1:nEq
        log10alpha(j) = log10(NND.vEta(j)) - mean(log10(kappa(j,:)));
    end
    alpha = 10.^log10alpha;
    for j = 1:nEq
        flag = true;
        if alpha(j)*A0 > 1.0
            flag = false;
        else
            if alpha(j)*A0 > rand()
                flag = false;
            end
        end
        if flag 
            NND.vClust(j) = 1;              % it belongs to a cluster 
        else                                % otherwise event j becomes a singleton event
            NND.vP_indx(j) = j;             % for a singleton event the index of a parent points to itself
            NND.vP_mag(j)  = NND.vCh_mag(j);
            NND.vDt(j)     = Inf;
            NND.vDist(j)   = Inf;
            NND.vClust(j)  = 0;             % it becomes a singleton event
        end
    end
end


function eta = nn_prox(eq_i,Ck,Par,sDistance)
%
    b = Par.b; df = Par.df; r_min = Par.r_min; r_max = Par.r_max;
    FMDi = 10.^(-b*eq_i(4));

    if strcmp(sDistance,'Radial')
        distance_j = @(lati,loni,latj,lonj) radialdist(lati,loni,latj,lonj);
    elseif strcmp(sDistance,'Euclidean')
        distance_j = @(yi,xi,yj,xj) sqrt((xi - xj)^2 + (yi - yj)^2);
    end

    nEqCk = length(Ck(:,1));
    mEta_j = Inf(1,nEqCk);
    for j = 1:nEqCk
        t_j = Ck(j,1) - eq_i(1); % j is a subsequent event
        %t_j = eq_i(1) - Ck(j,1); % j is a subsequent event
        r_j = distance_j(eq_i(2),eq_i(3),Ck(j,2),Ck(j,3));
        if t_j > 0.0 && r_j > r_min && r_j <= r_max % excludes events that are either occurred at the same place or time
            mEta_j(j) = t_j*(r_j)^df*FMDi;
        end
    end
    eta = min(mEta_j);
end
