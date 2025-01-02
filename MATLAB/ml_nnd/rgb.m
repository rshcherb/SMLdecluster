function col = rgb(sCol)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    if strcmp(sCol,'red')
        col = [1 0 0];
    elseif strcmp(sCol,'green')
        col = [0 1 0];
    elseif strcmp(sCol,'blue')
        col = [0 0 1];
    elseif strcmp(sCol,'black')
        col = [0 0 0];
    elseif strcmp(sCol,'gray')
        col = [0.5 0.5 0.5];
    end
end