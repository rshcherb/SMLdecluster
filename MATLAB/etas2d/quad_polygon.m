function [fUintR, DT, vXp, vYp] = quad_polygon(vX,vY,mU,vReg_poly,varargin)
%
%   Computes the quadrature of function mU over a polygon vReg_poly
%
%   vX, vY - linear coordinate vectors over which mU is defined
%   mU     - the function to integrate given as a matrix of size [length(vY) length(vX)]
%   vReg_poly - domain of integration as a polygon
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 15 May 2021
%   ...
%   version: 1.0.0, 15 May 2021
%
    sMethod = 'ZeroOutsidePolygon';
    for k = 1:length(varargin)
        if strcmp('Method',varargin{k})
            sMethod = varargin{k+1};
        end
    end

    [mX, mY] = meshgrid(vX,vY);
    inR = inpolygon(mX(:),mY(:),vReg_poly(:,1),vReg_poly(:,2)); % logical indeces of XY-coordinates inside the region R
    vXp = mX(inR);                                        % extract the grid points inside the polygon
    vYp = mY(inR);
    if strcmp(sMethod,'Triangulation')
        % using a triangulation
        DT = delaunay(vXp,vYp);                           % create a Delauney triangulation of the regular meshgrid inside the region R

        fUintR = 0;
        for nt = 1:length(DT(:,1))
            p1 = [vXp(DT(nt,1)), vYp(DT(nt,1))];          % three points of each triangle that form the triangulation
            p2 = [vXp(DT(nt,2)), vYp(DT(nt,2))];
            p3 = [vXp(DT(nt,3)), vYp(DT(nt,3))];
            Qnp = quadtriangle(8,'Type','nonproduct','Domain',[p1; p2; p3]); % for each triangle compute the quadrature
            vU = interp2(vX,vY,mU,Qnp.Points(:,1),Qnp.Points(:,2));
            fUintR = fUintR + dot(Qnp.Weights,vU);
        end
    elseif strcmp(sMethod,'ZeroOutsidePolygon')
        % using trapz() but setting the values of mU() outside region R to zero
        mU_inR = zeros(length(vX),length(vY));
        mU_inR(inR) = mU(inR);
        fUintR  = trapz(vY,trapz(vX,mU_inR,2));        % integral of the intensity over the whole area
        DT = [];
    end
end
