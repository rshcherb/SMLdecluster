function [fMstar] = gr_bath(fMc,vGRpar,vGRparErr)
%
% compute the inferred largest aftershock
% Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, November 9, 2014
%
    % computing the m^star of the modified Bath law
    fAval = vGRpar(2);
    fBest = vGRpar(1);
    fMstar = fAval/fBest;
end
