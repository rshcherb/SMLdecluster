function dt = mol_inverse(c,p)
%
%   Computes the time jump using the Omori kernel acording to Zhuang et al. (2005)
%   normalization: (p-1)/c*(1 + t/c)^(-p)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 April 2021
%   ...
%   version: 1.0.0, 28 April 2021
%
    U = 0.0;
    while U == 0.0
        U = rand();
    end
    
    dt = c*(U^(1/(1-p)) - 1);
end
