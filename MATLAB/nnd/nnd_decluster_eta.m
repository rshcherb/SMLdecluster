function NND = nnd_decluster_eta(NND,vEtaThresh,dt_max,dr_max)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.2.0, 3 December 2024
%
    nEq = size(NND.vEta,1);
    nthres = length(vEtaThresh);
    for j = 1:nEq
        for nth = 1:nthres
            %if NND.vEta(j) <= vEta_thresh(nth)
            if NND.vEta(j) <= vEtaThresh(nth) && NND.vDt(j) <= dt_max && NND.vDist(j) <= dr_max
                NND.vClust(j) = nth;            % it belongs to mode cluster n assuming there are several modes
            else                                % otherwise event j becomes a singleton event
                NND.vP_indx(j) = j;             % for a singleton event the index of a parent points to itself
                NND.vP_mag(j)  = NND.vCh_mag(j);
                NND.vDt(j)     = Inf;
                NND.vDist(j)   = Inf;
                NND.vClust(j)  = 0;             % it becomes a singleton event
            end
        end
    end
    NND.vEtaThresh = vEtaThresh;
    NND.vLog10EtaThresh = log10(vEtaThresh);
end

