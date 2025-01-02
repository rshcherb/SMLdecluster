function [V, C, A, vCircP] = voronoi_cell_area(vPoints,vReg)
%
%   Computes the Voronoi tesselation given vPoints inside a given region vReg. 
%   It also returns the area size of each voronoi cell
%
%   vPoints - vector of given points in space: [x, y]
%   vReg    - vector of region of interest: [x, y]
%   V       - list of all the vertices
%   C       - connectivity array: each row contains the list of vertex indeces forming a cell inside region vReg
%   A       - area size of each cell
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 25 December 2021
%   ...
%   version: 1.0.0, 27 December 2021
%
    nphi    = 20;                  % number of dummy points ouside the region vReg
    nP      = length(vPoints(:,1)); % number of original points
    vCircP  = zeros(nphi,2);       % dummy points on a circle sorrounding the region vReg
    xmin    = min(vReg(:,1));
    xmax    = max(vReg(:,1));
    ymin    = min(vReg(:,2));
    ymax    = max(vReg(:,2));
    fR = 2*max([xmax-xmin,ymax-ymin]); % radius of the circle where dummy points are placed
    x0 = 0.5*(xmax + xmin);        % x coordinates of the centre
    y0 = 0.5*(ymax + ymin);
    t = linspace(0,2*pi,nphi+1);
    t = t(1:end-1);
    t = t + 0.1*(0.5 - rand(1,nphi)); % add some noise to the locations of the points on a circle
    vCircP(:,1) = x0 + fR*cos(t); % x coordinates of a circle
    vCircP(:,2) = y0 + fR*sin(t); % y coordinates of a circle
    vPoints = [vPoints; vCircP];
    
    %vPoints = unique(vPoints,'rows','stable');  % eliminate non-unique points
    [V, C] = voronoin([vPoints(:,1),vPoints(:,2)],{'Qs'});

    poly_reg = polyshape(vReg(:,1),vReg(:,2));  % polygon object for the region
    nC = length(C);
    A  = zeros(nC,1);
    bI = false(nC,1);
    for i = 1:nC
        v1 = V(C{i},1); 
        v2 = V(C{i},2);
        if all(v1 ~= Inf) % only consider cells with verteces not having Inf values
            inR = inpolygon(v1,v2,vReg(:,1),vReg(:,2)); % find the vertices that are inside or outside vReg
            if all(inR)   % the whole cell is inside vReg
                A(i) = polyarea(v1,v2);
                bI(i) = true;
            else
                poly_i = polyshape(v1,v2);
                poly_int = intersect(poly_i,poly_reg);
                if poly_int.NumRegions == 1
                    A(i) = polyarea(poly_int.Vertices(:,1),poly_int.Vertices(:,2));
                    bI(i) = true;
                end
            end
        end
    end
    %V = V(bI);
    %C = C(bI);
    %A = A(bI);
    C = C(1:nP);
    A = A(1:nP);
end
