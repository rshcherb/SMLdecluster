%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 17 September 2024
%   ...
%   version 1.3.0, 19 November 2024
%
NNDpar.MixtModel = 'GMM';
NNDpar.ThreshMethod = 'minsaddle';
NNDpar.sEtaDistrbFit = 'Weibull'; % the model to fit to each mode of \eta distribution: 'Weibull', 'Normal'
if strcmp(Model.sRegion,'southcal')
    NNDpar.b = 1.04;
    NNDpar.df = 1.6;
    NNDpar.EtaLim = [10^(-10.5), 10^(1.5)];
    NNDpar.EtaYLim = [0, 0.7];
    NNDpar.TLim = [10^(-8.5), 10^(3.0)];
    NNDpar.RLim = [10^(-6.5), 10^(3.5)];
    %
    % NNDpar.b = 0.0;
    % NNDpar.EtaLim = [10^(-7), 10^(5)];
    % NNDpar.EtaYLim = [0, 0.8];
    % NNDpar.TLim = [10^(-6), 10^(5.0)];
    % NNDpar.RLim = [10^(-4), 10^(5.5)];
elseif strcmp(Model.sRegion,'italy')
    NNDpar.b = 1.1;
    NNDpar.df = 1.6;
    NNDpar.EtaLim = [10^(-10), 10^(2.5)];
    NNDpar.EtaYLim = [0, 0.65];
    NNDpar.TLim = [10^(-8), 10^(3.0)];
    NNDpar.RLim = [10^(-6), 10^(4.0)];
end
NNDpar.q = 0.5;
NNDpar.r_min = 0.0;
NNDpar.r_max = 1e10;
NNDpar.dt_max = 1e9;      % maximum time in days betweet parent-daughter pairs
NNDpar.dr_max = 1e10;     % maximum distance in m betweet parent-daughter pairs
NNDpar.nGMMmax = 2;
NNDpar.nNode_min = 50;    % the minimum number of nodes above which to plot
NNDpar.nClst_largest = 9; % the number of the largest clusters to plot
% duplicate in Model
Model.sEtaDistrbFit = NNDpar.sEtaDistrbFit; % the model to fit to each mode of \eta distribution: 'Weibull', 'Normal'
Model.EtaLim = NNDpar.EtaLim;
Model.EtaYLim = NNDpar.EtaYLim;
Model.TLim = NNDpar.TLim;
Model.RLim = NNDpar.RLim;
Model.nGMMmax = NNDpar.nGMMmax;
Model.nNode_min = NNDpar.nNode_min;    % the minimum number of nodes above which to plot
Model.nClst_largest = NNDpar.nClst_largest; % the number of the largest clusters to plot

Model.NNDpar = NNDpar;

