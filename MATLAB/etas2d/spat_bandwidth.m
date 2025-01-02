function vHj = spat_bandwidth(vP,nP,h_min,ns)
%
%   Computes the variable bandwidth for each event (Zhuang et al., JGR 2005, p.4, Eq. (21); Zhuang, EPS, 2011, Eq.(10))
%
%   vP    - [x, y] locations of the spatial points
%   nP    - number of points
%   h_min - the minumum bandwidth
%   ns    - the namber of nearest neighbours events to consider: integer in the typical range 3-10
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 19 May 2021
%   ...
%   version: 1.1.0, 29 December 2021
%
    %nP = length(vP(:,1));
    vHj = zeros(nP,1);
    for j = 1:nP % go over all the points
        dist_ij = zeros(nP,1);
        for i = 1:nP
            dist_ij(i) = sqrt( (vP(j,1) - vP(i,1)).^2 + (vP(j,2) - vP(i,2)).^2 ); % compute the distances between jth point and all the previous points
        end
        dist_ij = dist_ij(dist_ij > 0);
        [dist_min_np, indx_np] = mink(dist_ij,ns);    % find np shortest distances
        hmax = max(dist_ij(indx_np));                 % find the maximum among np distances. This defines the np_th nearest neighbour of the event j
        vHj(j) = max([h_min,hmax]);                   % the bandwidth is equal to either hj or the lower cutoff h_min
    end
end

