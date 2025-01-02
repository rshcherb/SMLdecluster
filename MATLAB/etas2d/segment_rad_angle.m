function [vRk, vRk2, vDk, vReg_ref, nVer_ref] = segment_rad_angle(vReg,inR,vEQxy,nNk,varargin)
%
%   Subdivides the given region, vReg, into segments originating from points inside/outside the region 
%
%   vReg     - the original region R bounding the target events: list of vertices in counterclockwise order vReb = [x, y]
%   inR      - the logical indeces of events in the catalogue that are in he region R
%   vEQxy    - the events locations: [x, y] or [lon, lat]
%   nNk      - the number of the initial radial segments
%
%   vRk      - the mean radius of each radial segment: size [nEQ, nVert_new]
%   vRk2     - for outside events this the distance to the second itercept of linae with the polygon R
%   vDk      - the angle of each radial segement in radians
%   vReg_new - a new refined list of vetices that include the original vReg vertices
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 7 May 2021
%   ...
%   version: 1.0.0, 11 June 2021
%
    sMethod = 'EqualDist'; % 'EqualAngle';
    bPlot   = false;
    for k = 1:length(varargin)
        if strcmp('Method',varargin{k})
            sMethod = varargin{k+1};
        end
    end
    
    nVer = length(vReg(:,1)) - 1; % the last vertex of vReg coincides with the first vertex
    if strcmp(sMethod,'EqualDist')
        % this is for events inside the target region R. The radial
        % segments are computed by connecting each event with the points on the boundary of the target region R 
        % create the list of new vertices along the polygon boundary
        [vReg_ref, nVer_ref] = poly_refine_vert(vReg,nNk); % nVer_ref may differ from the initial value nNk
        %disp([nNk,nVer_ref])
        nNk  = nVer_ref;
        nEq  = length(vEQxy(:,1));
        vRk  = zeros(nEq,nNk);
        vRk2 = zeros(nEq,nNk); % only for the the events outside the region R
        vDk  = zeros(nEq,nNk);
        for n = 1:nEq
            if inR(n) % only for events inside the target region R
                for k = 1:nVer_ref
                    a = sqrt( (vEQxy(n,1)    - vReg_ref(k,1))^2   + (vEQxy(n,2)    - vReg_ref(k,2))^2 );
                    b = sqrt( (vEQxy(n,1)    - vReg_ref(k+1,1))^2 + (vEQxy(n,2)    - vReg_ref(k+1,2))^2 );
                    c = sqrt( (vReg_ref(k,1) - vReg_ref(k+1,1))^2 + (vReg_ref(k,2) - vReg_ref(k+1,2))^2 );
                    vDk(n,k) = acos((a^2 + b^2 - c^2)/(2*a*b)); % cosine rule for computing the angle between sides a and b with the oposite side c
                    vRk(n,k) = 0.5*(a + b);                     % the radius is simply the average between sides a and b
                end
            else      % for events outside R
                % find two vertices from the original vReg connected to the event n that have the largest angle between the corresponding lines
                anglemax = 0;
                fRmax = 0;
                for i = 1:nVer-1
                    a = sqrt( (vEQxy(n,1) - vReg(i,1))^2 + (vEQxy(n,2) - vReg(i,2))^2 );
                    if a > fRmax
                        fRmax = a;
                    end
                    for j = i+1:nVer  % going counterclockwise from vertex i
                        vi = [vReg(i,1) - vEQxy(n,1), vReg(i,2) - vEQxy(n,2)]; % vector connecting event location with vertex i
                        vj = [vReg(j,1) - vEQxy(n,1), vReg(j,2) - vEQxy(n,2)]; % vector connecting event location with vertex j
                        angle = atan2(vi(1)*vj(2)-vi(2)*vj(1), vi(1)*vj(1)+vi(2)*vj(2)); % angle between vectors vi and vj
                        %disp(angle)
                        if abs(angle) > anglemax
                            if angle >= 0
                                v1x = vReg(i,1);
                                v1y = vReg(i,2);
                                v2x = vReg(j,1);
                                v2y = vReg(j,2);
                            else
                                v1x = vReg(j,1);
                                v1y = vReg(j,2);
                                v2x = vReg(i,1);
                                v2y = vReg(i,2);
                            end
                            anglemax = abs(angle);
                        end
                    end
                end
                fRmax = 2*fRmax;
%                 disp([v1x,v1y])
%                 disp([v2x,v2y])
                vPmax      = zeros(nNk,2);
                vPmax(1,:) = [v1x,v1y];
                vPmin      = zeros(nNk,2);
                vPmin(1,:) = [v1x,v1y];
                v0   = [1, 0];
                v1   = [v1x - vEQxy(n,1), v1y - vEQxy(n,2)];
                phi0 = atan2(v0(1)*v1(2)-v0(2)*v1(1), v0(1)*v1(1)+v0(2)*v1(2)); % angle between the vector v1 and the unit vector v0 along x-axis
                %disp(rad2deg(dphi))
                phi  = linspace(phi0,anglemax+phi0,nNk); % phi is the angles netween two end vertices separated by anglemax
                %disp(rad2deg(phi))
                for k = 2:nNk-1
                    % find the sides of a tringle formed by the event outside the region vReg, the point (v1x,v1y) and the
                    % line rotated counterclockwise by angle phi(k) and distanse fRmax away from the event location
                    
                    xk = vEQxy(n,1) + fRmax*cos(phi(k));
                    yk = vEQxy(n,2) + fRmax*sin(phi(k));
                    %disp([xk,yk])
                    
                    % this computes the intercept points of the line with the polygon vReg
                    [xi,yi] = polyxpoly([vEQxy(n,1) xk],[vEQxy(n,2) yk],vReg(:,1),vReg(:,2));
                    %disp(k)
%                     disp([vEQxy(n,1),vEQxy(n,2)])
                    %disp([xi,yi])
                    r1 = sqrt( (vEQxy(n,1) - xi(1))^2 + (vEQxy(n,2) - yi(1))^2 ); % the distance between the event and the vertex
                    r2 = sqrt( (vEQxy(n,1) - xi(2))^2 + (vEQxy(n,2) - yi(2))^2 ); % the distance between the event and the vertex
                    if r1 < r2
                        r_min = r1;
                        r_max = r2;
                        vPmin(k,:) = [xi(1),yi(1)];
                        vPmax(k,:) = [xi(2),yi(2)];
                    else
                        r_min = r2;
                        r_max = r1;
                        vPmin(k,:) = [xi(2),yi(2)];
                        vPmax(k,:) = [xi(1),yi(1)];
                    end
                    vDk(n,k-1)  = phi(k) - phi(k-1);               % the angle between sides b and c with the oposite side a
                    vRk(n,k-1)  = 0.5*(r_min + vRk(n,k-1));        % the radius is simply the average between sides b and c
                    vRk2(n,k-1) = 0.5*(r_max + vRk2(n,k-1));       % the radius is simply the average between sides b and c
                end
                vPmin(nNk,:) = [v2x,v2y];
                vPmax(nNk,:) = [v2x,v2y];
                
                if bPlot
                    figure('Name','Polygon segments','Position',[550 400 800 600]);

                    pgon = polyshape(vReg);
                    vPj = ones(size(vPmax)).*vEQxy(n,:);
                    mX = [vPj(:,1)'; vPmax(:,1)']; % columns of x for the point Pj and points on the boundary of the polygon
                    mY = [vPj(:,2)'; vPmax(:,2)']; % columns of y for the point Pj and points on the boundary of the polygon

                    plot(vReg(:,1),vReg(:,2),'o');
                    hold on
                    plot(vPmax(:,1),vPmax(:,2),'ks');    % points along the boundary of the polygon
                    plot(vPmin(:,1),vPmin(:,2),'kx');    % points along the boundary of the polygon
                    plot(pgon);
                    plot(vEQxy(n,1),vEQxy(n,2),'d');     % point Pj outside the polygon
                    plot(vReg(1,1),vReg(1,2),'rd');      % the first point in vReg
                    text(vReg(1,1),vReg(1,2),'1');      % the first point in vReg
                    plot(mX,mY);
                    %plot(v1x,v1y,'rx',v2x,v2y,'bx');
                    text(v1x,v1y,'v1')
                    text(v2x,v2y,'v2')

                    ylabel('lat');
                    xlabel('lon');
                    hold off
                end
            end
        end
    elseif strcmp(sMethod,'EqualAngle')

    end
end
