function [h, p, tstat] = ks_test(vTimes,options)
%
%   The KS test asseses the deviation of the cumulative number of
%   events from what is expected in a stationary Poisson process with the
%   same number of total events N (Zaliapin and Ben-Zion, JGR, 2020, p12)
%   vTimes - event occurrence times
%   h - logical: 'false' or 0 if the test fails to reject the null hypothesis
%   p - p-value of the test
%   tstat - test statistic
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 27 October 2024
%   ...
%   version 1.1.0, 30 October 2024
%
    arguments
        vTimes double
        options.Alpha double = 0.05  % default 5% significance level
        options.Display char = 'off'    % whether to disply the result of the test
    end
    tmax = max(vTimes);           % the time of the first event may occur for t > 0
    tmin = min(vTimes);           % the time of the last event 

    u_i = (vTimes - tmin)/(tmax - tmin);
    test_cdf = makedist("Uniform");
    [h, p, tstat, cv] = kstest(u_i,'CDF',test_cdf,'Alpha',options.Alpha);

    % if p < options.Alpha
    %     h = true;       % the hypothesis that the event times come from the Poisson process is rejected
    % else
    %     h = false;      % test fails to reject the hypothesis at the given significance level OptArgs.Alpha
    % end
    if strcmp(options.Display,'on')
        disp(['p-value: ',num2str(p),'; test statistic: ',num2str(tstat)])
        if h
            disp('The KS test has rejected the null hypothesis (Poisson distributed)')
        else
            disp('The KS test has failed to reject the null hypothesis (Poisson distributed)')
        end
    end
end

