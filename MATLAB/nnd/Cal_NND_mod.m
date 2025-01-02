function [NND,EQEvents] = Cal_NND_mod(mCat,vPar)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Input
%   mCat [time lan lon mag depth timeyear]
    EQEvents = table;
    EQEvents.OriginTime = mCat(:,1);
    EQEvents.Latitude = mCat(:,2);
    EQEvents.Longitude = mCat(:,3);
    EQEvents.Mag = mCat(:,4);
    %EQEvents.Depth = mCat(:,5);
    EQEvents = sortrows(EQEvents,1);
    
     % Set constants as per Zaliapin et al., 2013
    w = 1;
    b    = vPar(1);
    df   = vPar(2);
    q    = vPar(3);
    nEq = length(EQEvents.Mag);
    NND = zeros(nEq,10);
    %disp(nEq)
    for ii = nEq:-1:2
        child  = EQEvents(ii,:); % Set child event
        parent = EQEvents(EQEvents.Mag >= child.Mag(1),:); % Set parent event
        t = (child.OriginTime(1) - parent.OriginTime);
        % If t is >= 0, make inf
        t(t <= 0) = inf;
        % Get spatial distance
        r = distance_hvrsn(child.Latitude,child.Longitude,parent.Latitude,parent.Longitude);
        
        r(r == 0) = inf;
        n = t.*r.^(df).*10.^(w*-b.*parent.Mag); % eta from equation 1
        [~,ind] = min(n);
        % Get minimum and save values
        T(ii-1,1) = t(ind)*10^(w*-q*b*parent.Mag(ind));
        R(ii-1,1) = r(ind)^(df)*10^(w*-(1-q)*b*parent.Mag(ind));
        N(ii-1,1) = T(ii-1,1)*R(ii-1,1);%n(ind);
        indp = find(ismember(EQEvents,parent(ind,:))==1);
        %disp([indp,ind])
        NND(ii-1,:) = [indp(1),ii,t(ind),r(ind),n(ind),T(ii-1,1),R(ii-1,1),N(ii-1,1),parent.Mag(ind),child.Mag(1)];
    end

    fid = fopen('nnd_values_mohammad.dat','w');
    %fprintf(fid,'parent ind \n');
    for n = 1:nEq-2
        fprintf(fid,'%4d %4d %13f %11f %11f %11f %11f %11f %11f\n',...
            NND(n,1),NND(n,2),NND(n,3),NND(n,4),NND(n,6),NND(n,7),NND(n,8),NND(n,9),NND(n,10));
    end
    fclose(fid);
end

function [dist_km] = distance_hvrsn(lat1, lon1, lat2, lon2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the distance between two points on the globe using the
% haversine formula.
%
% Input:
% lat1   - decimal latitude of the first point
% lon1   - decimal longitude of the first point
% depth1 - Depth of the first point
% lat2   - decimal latitude of the second point
% lon2   - decimal longitude of the second point
% 
% Output
% dist_degree - distance in degrees 2d
% dist_km - distance in kilometers 3d
% hrs_dist_km distance in kilometers 2d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R =  6371000;    % Earth's radius (m)
d2r = pi/180;    % deg = rad*180/pi
% Convert degree latitude and longitude to radians
si1 = d2r*(lat1);
si2 = d2r*(lat2);
% Convert latitude and longitude differences to radians 
si_del = d2r*(lat2 - lat1);
lam_del = d2r*(lon2 - lon1);
% Haversine forumlation for distance
a = (sin(si_del/2).^2) + (cos(si1).*cos(si2).*(sin(lam_del/2)).^2);
c = 2 * atan2(sqrt(a), sqrt(1-a));
%Convert distance radians back to degrees
dist_km = (R*c)/1000;
%dist_degree = c * 180/pi;
end

